
#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 2) stop("Incorrect number of arguments, usage> Rscript index.r SIZE OUTPUT"))
index=args[1]
opt_output = args[2]
a<-sample(1:407, size = index)
write.table(a, opt_output,row.names=F,col.names=F) 
