#!/zs32/home/frwang/software/R-3.4.0/bin/Rscript
library(qvalue)
for(i in c(1:22)){
    d = read.table(paste("sqtl.nopermute.chr",i,sep=""), head=F, sep=" ", stringsAsFactors=F)
    d$qvalue = qvalue(d$V12)$qvalues
    write.table(d[which(!is.na(d$qvalue)),], paste("sqtl.nopermute.chr",i,".qvalue",sep=""), quote=F, row.names=F, col.names=T, sep="\t")
}
