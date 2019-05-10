#Load package
suppressMessages(library(qvalue))
suppressMessages(library(limma))
suppressMessages(library(DESeq2))
suppressMessages(library(WGCNA))
suppressMessages(library(preprocessCore))
suppressMessages(library(peer))

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("Incorrect number of arguments, usage> Rscript RNA_processing.R INPUT PEER_factor OUTPUT"))
opt_input  = args[1];
opt_factor = as.numeric(args[2]);
opt_outputdir = args[3];

count<-read.table(opt_input,row.names=1,stringsAsFactors=F)
colnames(count)<-count[1,]
count2<-read.table(opt_input,row.names=1,stringsAsFactors=F,head=T)
colnames(count2)<-colnames(count)
log2cpm<-voom(count2)$E
genes_to_keep = apply(log2cpm>=0,1,sum) >= round(0.25 * ncol(log2cpm))
annotation = read.csv("/zs32/data-analysis/liucy_group/shareData/Chinese_brain/RNA-seq/expression/annotation.gene.gencodeV19.csv") 
log2cpm.fgene = log2cpm[genes_to_keep,] 
normadj <- (0.5+0.5*bicor(log2cpm.fgene, use='pairwise.complete.obs'))^2
#normadj <- (0.5+0.5*cor(log2cpm.fgene))^2
netsummary <- fundamentalNetworkConcepts(normadj);
k <- netsummary$Connectivity; K <- k/max(k); Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
sdout <- 3
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(log2cpm.fgene)[outliers]); print(table(outliers))
log2cpm.fgene.fsample <- log2cpm.fgene[,!outliers]
pdf(paste(opt_outputdir,"/log2cpm.fgene.zscore.pdf",sep=""))
#plot(Z.K, col = as.numeric(colData(dds)$Region), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
plot(Z.K, col = "black", pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-sdout, lty=2)
dev.off()
log2cpm.fgene.fsample.qn<-normalize.quantiles(log2cpm.fgene.fsample,copy=T)  # Quantile normalization across columns
rownames(log2cpm.fgene.fsample.qn)<-rownames(log2cpm.fgene.fsample)
colnames(log2cpm.fgene.fsample.qn)<-colnames(log2cpm.fgene.fsample)
out=paste(opt_outputdir, "/log2cpm.fgene.fsample.qn", sep="")
write.table(log2cpm.fgene.fsample.qn, out,sep="\t", row.names=T, quote=F,col.names=NA)
expr = t(as.matrix(log2cpm.fgene.fsample.qn))  # N rows and G columns, where N is the number of samples, and G is the number of genes. No column or row names.
dim(expr)
factorList=list()
residList=list()
for(nFactor in 1:1*opt_factor){
    model = PEER()  # create the model object
    PEER_setPhenoMean(model,expr)  # set the observed data
    #PEER_setNk(model,30)  # Set the hidden confounders. 
    PEER_setNk(model,nFactor) # gradient number of factors
    PEER_getNk(model)
    PEER_setAdd_mean(model, TRUE)  # include an additional factor (covariate) to account for the mean expression
    ## PEER_setCovariates(model, as.matrix(meta))  # adding covariates has no effect on the model?
    #PEER_setNmax_iterations(model, 100)  # If the model is not converged after the default 1,000 iterations, and the variance of the residuals keeps decreasing, choose a higher value of iterations, e.g., 10,000.
    PEER_update(model)  # perform the inference
    factors = PEER_getX(model)  # inferred confounders
    weights = PEER_getW(model)  # their weights
    precision = PEER_getAlpha(model)     # precision (inverse variance) of the weights
    residuals = PEER_getResiduals(model) # the residual dataset
    #plot(precision)
    #PEER_plotModel(model)
    Variance<-c(); for(i in 2:(opt_factor+1)){Variance<-c(Variance,var(factors[,i]))}; Variance<-sort(Variance, decreasing=TRUE); Variance<-100*Variance/sum(Variance)
    pdf(paste(opt_outputdir,"/variance-factor.pdf",sep=""),width=5,height=5)
    plot(Variance,type="o",pch=19,col="red",xlab="Factors")
    dev.off()
    factors1<-t(factors)[-1,]
    colnames(factors1)<-colnames(log2cpm.fgene.fsample.qn)
    rownames(factors1)<-paste("Factor",1:nrow(factors1),sep="")
    factorList[[nFactor]]<-factors1
    #factors2<-rbind(t(meta),factors1)
    #factors2<-cbind(meta,t(factors1))
    #write.table(t(factors2), file=paste("log2cpm.fgene.fsample.qn.realign.peerCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
    write.table(factors1, file=paste(opt_outputdir,"/expr",opt_factor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
}
