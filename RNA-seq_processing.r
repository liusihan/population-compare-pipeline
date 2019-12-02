###########################################
##-------- Read meta information --------##
##-------- Cluster library batch --------##
###########################################

## load meta information
meta<-read.table("meta.txt",sep="\t",head=T)
meta$LibraryBatch<-as.character(meta$LibraryBatch)

## LibraryBatch cluster by hclust
LibraryBatch<-matrix(as.numeric(as.Date(meta$LibraryBatch)), ncol=1)
rownames(LibraryBatch)<-meta$LibraryBatch
LibraryBatch<-unique(LibraryBatch)
LibraryBatch.dist<-dist(LibraryBatch,method="euclidean")
LibraryBatch.hclust<-hclust(LibraryBatch.dist, method="median")
LibraryBatch.id<-cutree(LibraryBatch.hclust,k=18)
pdf("hclust.LibraryBatch.pdf",width=12,height=8)
plot(LibraryBatch.hclust, main="Cluster Dendrogram", xlab="Labrary batch", ylab="Euclidean distance")
rect.hclust(LibraryBatch.hclust,k=18)
dev.off()

## RNAIsolationBatch cluster by hclust
RNAIsolationBatch<-matrix(as.numeric(as.Date(meta$RNAIsolationBatch)), ncol=1)
rownames(RNAIsolationBatch)<-meta$RNAIsolationBatch
RNAIsolationBatch<-unique(RNAIsolationBatch)
RNAIsolationBatch.dist<-dist(RNAIsolationBatch,method="euclidean")
RNAIsolationBatch.hclust<-hclust(RNAIsolationBatch.dist, method="median")
RNAIsolationBatch.id<-cutree(RNAIsolationBatch.hclust,k=14)
pdf("hclust.RNAIsolationBatch.pdf",width=12,height=8)
plot(RNAIsolationBatch.hclust, main="Cluster Dendrogram", xlab="Labrary batch", ylab="Euclidean distance")
rect.hclust(RNAIsolationBatch.hclust,k=14)
dev.off()

## output formated meta information
meta$libraryBatchClustered<-rep("",nrow(meta))
for(i in 1:nrow(meta)){
    meta$libraryBatchClustered[i]<-paste("B",LibraryBatch.id[meta$LibraryBatch[i]],sep="")
}
meta$RNAIsolationBatchClustered<-rep("",nrow(meta))
for(i in 1:nrow(meta)){
    meta$RNAIsolationBatchClustered[i]<-paste("B",RNAIsolationBatch.id[meta$RNAIsolationBatch[i]],sep="")
}
rownames(meta)<-meta$BID
write.table(meta, file="meta2.txt", sep="\t", row.names=F, col.names=T, quote=F)

meta$RIN_square <- meta$RIN^2
meta$ageDeath_square <- meta$ageDeath^2

###########################################
##--- Read RNA quantification data   ----##
##--- log2(CPM) normalization        ----##
###########################################

#read raw count data and annotation file
count<-read.table("raw.count",head=T,row.names=1)
annotation = read.csv("annotation.gene.gencodeV19.csv")

# assess for Sample Swaps by XIST and Y chromosome
mds = cmdscale(dist(t(log2(0.001 + count[annotation$chr=="chrY",]))))
xist = log2(0.001+count["ENSG00000229807",])

col.blue = rgb(t(col2rgb("blue")),alpha=50,maxColorValue = 255); col.pink=rgb(t(col2rgb("pink")),alpha=50,maxColorValue = 255)
sex_col = rep(col.blue, times=nrow(mds)); sex_col[meta$Sex=="F"] =col.pink
plot(xist, mds[,1], col=sex_col,pch=19)

tree = hclust(dist(cbind(xist,mds[,1])),"average")
sex_pred= factor(gsub(2, "F", gsub("1","M", cutree(tree,k=2))))
plot(xist, mds[,1], col=(sex_pred),pch=19)

discordant = (meta$sex=="M" & sex_pred=="F") | (meta$sex=="F" & sex_pred=="M")

#remove sex mismatch sample
count2<-count[,!discordant]


#log2(CPM) normalization
library(limma)
log2cpm<-voom(count2)$E

write.table(log2cpm,file="log2cpm",sep="\t",row.names=T,col.names=T,quote=F)

###########################################
##---  filter                        ----##
###########################################

## keep genes with at least 1 CPM in at least 25% of the individuals
library(DESeq2)
genes_to_keep = apply(log2cpm>=0,1,sum) >= round(0.25 * ncol(log2cpm))
log2cpm.fgene = log2cpm[genes_to_keep,] 
annot<-annotation[annotation$gene_id %in% rownames(log2cpm.fgene),]
brainExpressedNoMT<-annot[annot$chr!="chrM"&annot$chr!="chrY"&annot$chr!="chrX",]
log2cpm.fgene = log2cpm[brainExpressedNoMT$gene_id,]
write.table(log2cpm.fgene, file="log2cpm.fgene", sep="\t", row.names=T, quote=F,col.names=T)

## pca
library(ggfortify)
library(scales)
pdf("log2cpm.fgene.pca.pdf")
pc<-prcomp(t(log2cpm.fgene))
#autoplot(pc)
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
eigsP1<-percent(eigs[1]/sum(eigs))
eigsP2<-percent(eigs[2]/sum(eigs))
plot(pc$x[,1],pc$x[,2],pch=19,col="grey",main="PCA",xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()

## Remove outlier samples
## Code modified from Michael Gandal's Github: https://github.com/mgandal/TSC_MIA_RNAseq/blob/master/code/step4a_Expression_Analysis.R
## Original paper: Network methods for describing sample relationships in genomic datasets: application to Huntington’s disease
library(WGCNA)
normadj <- (0.5+0.5*bicor(log2cpm.fgene, use='pairwise.complete.obs'))^2
#normadj <- (0.5+0.5*cor(log2cpm.fgene))^2
netsummary <- fundamentalNetworkConcepts(normadj);
k <- netsummary$Connectivity; K <- k/max(k); Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
sdout <- 3
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(log2cpm.fgene)[outliers]); print(table(outliers))
color<-rownames(pc$x) %in% colnames(log2cpm.fgene)[outliers]
color[which(color==FALSE)]<-"grey"
color[which(color==TRUE)]<-"red"
pdf("log2cpm.fgene.zscore.pdf")
#plot(Z.K, col = as.numeric(colData(dds)$Region), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
plot(Z.K, col = "black", pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-sdout, lty=2)
dev.off()
pdf("log2cpm.fgene.pca.markOutlier.pdf")
plot(pc$x[,1],pc$x[,2],pch=19,col=color,main=paste("PCA (marked Z-score < -",sdout,")",sep=""),xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()
log2cpm.fgene.fsample <- log2cpm.fgene[,!outliers]
write.table(log2cpm.fgene.fsample, file="log2cpm.fgene.fsample", sep="\t", row.names=T, quote=F,col.names=T)


###########################################
##--- Quantile Normalization       ------##
###########################################

## QN
library(preprocessCore)
log2cpm.fgene.fsample.qn<-normalize.quantiles(log2cpm.fgene.fsample,copy=T)  # Quantile normalization across columns
rownames(log2cpm.fgene.fsample.qn)<-rownames(log2cpm.fgene.fsample)
colnames(log2cpm.fgene.fsample.qn)<-colnames(log2cpm.fgene.fsample)
write.table(combat_edata, file="log2cpm.fgene.fsample.qn.combat", sep="\t", row.names=T, quote=F,col.names=NA)


###
#Batch Clustering
pdf("Batch_Clustering.pdf",width=15,height=10)
PC = prcomp(na.omit(t(scale(t(log2cpm.fgene.fsample.qn),scale=F))),scale=F)
PC = PC$rotation[,1:5]
library(WGCNA)
batch = model.matrix(~0+LibraryBatch+RNAseqBatch+RNAIsolationBatch, data=meta)
idx = match(rownames(batch), rownames(meta))
tree=hclust(dist((batch)),method="ward.D2")
colors=with(meta[idx,], cbind(labels2colors(LibraryBatch), 
                                 labels2colors(RNAseqBatch), 
                                 labels2colors(RNAIsolationBatch)))
                                 
colors = cbind(colors, numbers2colors(PC[idx,]))            
plotDendroAndColors(tree,colors, groupLabels = c("LibraryBatch", "RNAseqBatch", "RNAIsolationBatch", paste0("PC",1:5)),cex.dendroLabels = .7,main="Batch Clustering")


cutsteps = seq(min(tree$height), 0.98*max(tree$height), by = (max(tree$height) - min(tree$height))/50)
r2out = numeric()
for(h in cutsteps) {
  c=as.factor(cutree(tree,h = h))
  s=0
  for(i in 1:5) s=s+summary(lm(PC[idx,1]~c))$adj.r.squared/5
  r2out=c(r2out,s)
}
plot(cutsteps, r2out)

dev.off()


## pca for quantile data
pdf("log2cpm.fgene.fsample.qn.pca.pdf")
pc<-prcomp(t(log2cpm.fgene.fsample.qn))
#autoplot(pc)
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
percVar<-eigs/sum(eigs)
eigsP1<-percent(eigs[1]/sum(eigs))
eigsP2<-percent(eigs[2]/sum(eigs))
plot(pc$x[,1],pc$x[,2],pch=19,col="grey",main="PCA (Quantile Normalized)",xlab=paste("PC1 (",eigsP1,")",sep=""),ylab=paste("PC2 (",eigsP2,")",sep=""))
dev.off()

###########################################
##---- Estimate hidden factors        ---##
###########################################

library(peer)
expr = t(as.matrix(log2cpm.fgene.fsample.qn))  # N rows and G columns, where N is the number of samples, and G is the number of genes. No column or row names.
dim(expr)
factorList=list()
residList=list()


for(nFactor in 1:1*20){
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
    Variance<-c(); for(i in 2:21){Variance<-c(Variance,var(factors[,i]))}; Variance<-sort(Variance, decreasing=TRUE); Variance<-100*Variance/sum(Variance)
    pdf("variance-factor.pdf",width=5,height=5)
    plot(Variance,type="o",pch=19,col="red",xlab="Factors")
    dev.off()
    factors1<-t(factors)[-1,]
    colnames(factors1)<-colnames(combat_edata)
    rownames(factors1)<-paste("Factor",1:nrow(factors1),sep="")
    factorList[[nFactor]]<-factors1
    #factors2<-rbind(t(meta),factors1)
    factors2<-cbind(meta,t(factors1))
    #write.table(t(factors2), file=paste("log2cpm.fgene.fsample.qn.realign.peerCovariates.factor",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
    write.table(factors1, file=paste("expr",nFactor,".xls",sep=""), sep="\t", row.names=T, quote=F,col.names=NA)
}

## test correlation of hidden factors and covariates
library(gplots)
rsqList<-list()
pvalList<-list()
for(nFactor in 1:1*20){
    factors1<-factorList[[nFactor]]
    rsq<-matrix(nrow=nFactor,ncol=ncol(meta))
    pval<-matrix(nrow=nFactor,ncol=ncol(meta))
    rownames(rsq)<-rownames(factors1)
    rownames(pval)<-rownames(factors1)
    colnames(rsq)<-colnames(meta)
    colnames(pval)<-colnames(meta)
    for(m in 1:nrow(factors1)){
        for(n in 1:ncol(meta)){
            colnames(meta)[n]
            mod<-model.matrix( ~ meta[,n])
            res<-summary(lm(factors1[m,] ~ mod))
            r2<-res$adj.r.squared
            p<-pf(res$fstatistic[1],res$fstatistic[2],res$fstatistic[3],lower.tail=FALSE)
            rsq[m,n]<-r2
            pval[m,n]<-p
        }
    }
    
    write.table(rsq,file=paste("cor.hidden-known.rSquare.factor",nFactor,".xls",sep=""),sep="\t",quote=F,col.names=NA)
    write.table(pval,file=paste("cor.hidden-known.pValue.factor",nFactor,".xls",sep=""),sep="\t",quote=F,col.names=NA)
    Lab.palette=colorRampPalette(c('white','red'),space="Lab")
    pdf(paste("cor.hidden-known.rSquare.factor",nFactor,".pdf",sep=""))
    heatmap.2(rsq,main="R-square",col=Lab.palette,density.info="none",trace="none",dendrogram="none",Colv=FALSE,Rowv=FALSE,scale="none",margins=c(12,6))
    dev.off()
    pdf(paste("cor.hidden-known.pValue.factor",nFactor,".pdf",sep=""))
    heatmap.2(-log(pval),main="-log(Pvalue)",col=Lab.palette,density.info="none",trace="none",dendrogram="none",Colv=FALSE,Rowv=FALSE,scale="none",margins=c(12,6))
    dev.off()
    
    rsqList[[nFactor]]<-rsq
    pvalList[[nFactor]]<-pval
}


##################################################
##- Bayesian Information Criterion (BIC) score -##
##- Akaike   Information Criterion (AIC) score -##
##################################################

colnames(log2cpm.fgene.fsample.qn)<-rownames(factors2)
d<-cbind(t(log2cpm.fgene.fsample.qn),factors2)
nGene<-nrow(log2cpm.fgene.fsample.qn)
nSample<-nrow(factors2)



## about two days for this step
covBIC<-c()
for(i in 1:nGene){
    capture.output(resBIC<-step(lm(d[,i] ~ d$CauseDeath+d$sex + d$Group+d$ageDeath +  d$RNAIsolationBatch + d$RIN + d$LibraryBatch + d$RNAseqBatch + d$ageDeath_square + d$RIN_square + d$Factor1 + d$Factor2 + d$Factor3 + d$Factor4 + d$Factor5 + d$Factor6 + d$Factor7 + d$Factor8 + d$Factor9 + d$Factor10 + d$Factor11 + d$Factor12 + d$Factor13 + d$Factor14 + d$Factor15 + d$Factor16 + d$Factor17 + d$Factor18 + d$Factor19 + d$Factor20 + d$seqPC1 + d$seqPC2 + d$seqPC3 + d$seqPC4 + d$seqPC5 + d$seqPC6 + d$seqPC7 + d$seqPC8 + d$seqPC9 + d$seqPC10+d$seqPC11+d$seqPC12+d$seqPC13+d$seqPC14), k=log(nSample), direction="both"), file=paste("BIC/log.BIC.",i,sep=""))
    covBIC<-c(covBIC, names(attr(resBIC$terms,"dataClasses"))[-1])
}

# covAIC<-c()
# for(i in 1:nGene){
#     capture.output(resAIC<-step(lm(d[,i] ~ d$BrainBank + d$Hemisphere + d$PMI + d$BrainWeight + d$YearAutopsy + d$Sex + d$Ethnicity + d$AgeDeath + d$Diagnosis + d$TissueState + d$RNAIsolationBatchClustered + d$RIN + d$ERCC_Added + d$libraryBatchClustered + d$RNAseqBatch + d$SequencingPlatform + d$PMI_square + d$BrainWeight_square + d$YearAutopsy_square + d$AgeDeath_square + d$RIN_square + d$Factor1 + d$Factor2 + d$Factor3 + d$Factor4 + d$Factor5 + d$Factor6 + d$Factor7 + d$Factor8 + d$Factor9 + d$Factor10 + d$Factor11 + d$Factor12 + d$Factor13 + d$Factor14 + d$Factor15 + d$Factor16 + d$Factor17 + d$Factor18 + d$Factor19 + d$Factor20 + d$Factor21 + d$Factor22 + d$Factor23 + d$Factor24 + d$Factor25 + d$Factor26 + d$Factor27 + d$Factor28 + d$Factor29 + d$Factor30), direction="both"), file=paste("log/log.AIC.",i,sep=""))
#     covAIC<-c(covAIC, names(attr(resAIC$terms,"dataClasses"))[-1])
# }

# Draw
phe<-colnames(factors2)
pdf("covariates.bic.pdf",width=7,height=3)
props <- table(covBIC)[paste("d$", phe, sep="")]/nGene
bp <- barplot(props,  xlab = "Covariates", ylab = "Proportion of genes decreased BIC", main="Select covariates by decreased BIC", ylim= c(0,1),col = c("blue"), las=2, cex.axis=0.5, cex.lab=0.5, cex.main=0.7, axisnames=FALSE)
#axis(1, at = bp, labels = effectsNames, xlab = "Covariates", cex.axis = 0.5, las=2)  # vertical x-axis
text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = phe, xpd = TRUE,cex=0.4)  # rotate 45 x-axis
text(bp, props, labels = round(props, 3), srt = 45, pos=3, cex = 0.3) # place numbers on top of bars 
dev.off()

# pdf("covariates.aic.pdf",width=7,height=3)
# props <- table(covAIC)[paste("d$", phe3, sep="")]/nGene
# bp <- barplot(props,  xlab = "Covariates", ylab = "Proportion of genes decreased AIC", main="Select covariates by decreased AIC", ylim= c(0,1),col = c("blue"), las=2, cex.axis=0.5, cex.lab=0.5, cex.main=0.7, axisnames=FALSE)
# #axis(1, at = bp, labels = effectsNames, xlab = "Covariates", cex.axis = 0.5, las=2)  # vertical x-axis
# text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = phe3, xpd = TRUE,cex=0.4)  # rotate 45 x-axis
# text(bp, props, labels = round(props, 3), srt = 45, pos=3, cex = 0.3) # place numbers on top of bars 
# dev.off()



save.image(file="RNA_processing.cpm.all.RData")
write.table(factors2,file="meta_with_peer.txt, sep="\t", row.names=T, quote=F,col.names=NA)
