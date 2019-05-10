#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 2) stop("Incorrect number of arguments, usage> Rscript Zscore.R INPUT OUTPUT"))
opt_input  = args[1];
opt_output  = args[2];

a<-read.table(opt_input,row.names=1)
a$Z<-scale(a$V2)[,1]
write.table(a,opt_output,col.names=F,row.names=T)
