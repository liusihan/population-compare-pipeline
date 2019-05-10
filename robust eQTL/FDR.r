library(qvalue)
for(i in c(1:22,"X")){
    d = read.table(paste("eqtl.nopermute.chr",i,sep=""), head=F, sep=" ", stringsAsFactors=F)
    d$qvalue = qvalue(d$V12)$qvalues
    write.table(d[which(!is.na(d$qvalue)),], paste("eqtl.nopermute.chr",i,".qvalue",sep=""), quote=F, row.names=F, col.names=T, sep="\t")
}
