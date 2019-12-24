
#!/zs32/home/frwang/software/R-3.4.0/bin/Rscript
library(qvalue)
library(preprocessCore)
library(peer)
#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("Incorrect number of arguments, usage> Rscript RNA_processing.R INPUT PEER_factor OUTPUT"))
opt_input  = args[1];
opt_factor = as.numeric(args[2]);
opt_outputdir = args[3];
qnpsi<-read.table(opt_input,row.names=1,head=T,stringsAsFactors=F)
expr = t(as.matrix(qnpsi))  # N rows and G columns, where N is the number of samples, and G is the number of genes. No column or row names.
dim(expr)
factorList=list()
residList=list()
for(nFactor in 1:1*opt_factor){
    model = PEER()  # create the model object
    PEER_setPhenoMean(model,expr)  # set the observed data
    #PEER_setNk(model,30)  # Set the hidden confounders. 
    PEER_setNk(model,nFactor) # gradient number of factors
    PEER_getNk(model)
    PEER_setAdd_mean(model, TRUE)  # include an additional factor (covariate) to acqnpsi for the mean expression
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
    colnames(factors1)<-colnames(qnpsi)
    rownames(factors1)<-paste("Factor",1:nrow(factors1),sep="")
    write.table(factors1, file=paste(opt_outputdir,"/expr",opt_factor,".txt",sep=""), sep="\t", row.names=T, quote=F,col.names=F)
}
