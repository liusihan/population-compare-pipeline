#Load package
suppressMessages(library(qvalue))

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 1) stop("Incorrect number of arguments, usage> Rscript pi1.R INPUT"))
opt_input  = args[1];

 #caucalate pi1
a<-read.table(opt_input)
library(qvalue)
pi0<-pi0est(a$V1,lambda=seq(0.1,0.9,0.1),method="smoother")
print(1-pi0$pi0)
