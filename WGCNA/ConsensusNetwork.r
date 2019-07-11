library(WGCNA)
multiExpr = vector(mode="list", length=2)
multiExpr[[1]] = list(data=as.data.frame(t(CAU_datExpr[genes,])))
multiExpr[[2]] = list(data=as.data.frame(t(CNS_datExpr[genes,])))
#multiExpr[[3]] = list(data=as.data.frame(t(AA_datExpr[genes,])))

bsize = 5000
nSets = 2
#nSets = 3
powers = c(seq(1,19,by=1),seq(20,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()

#Pick Soft Threshold
if(TRUE) {  
  pdf("./unsigned_softThresh.pdf")
  for(n in 1:2) {
    multiExpr[[n]]$softThresh = pickSoftThreshold(data= multiExpr[[n]]$data, networkType = "unsigned", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)
  
  sft = multiExpr[[n]]$softThresh
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n", main=names(multiExpr)[n])
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2",main=names(multiExpr)[n])
  abline(h=0.8, col="black")
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
  }
  dev.off()
}


softPower=5; cquant=0.5

net = blockwiseConsensusModules(multiExpr, power = 9, maxBlockSize=50000, minModuleSize = 50, deepSplit = 2,pamRespectsDendro = FALSE,networkType="signed",mergeCutHeight = 0.3, numericLabels = TRUE,minKMEtoStay = 0,saveTOMs = TRUE, verbose = 5, pamStage=FALSE ,saveConsensusTOMs = TRUE,consensusTOMFileNames = "Signed_WGCNA_consensusTOM.RData")

								
consMEs = net$multiMEs;
moduleLabels = net$colors;
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]];
pdf("Unsigned gene dendrogram and module colors.pdf")
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Unsigned Consensus gene dendrogram and module colors")
dev.off()

save(net,consMEs, moduleColors, moduleLabels, consTree, file = "Signed_Consensus-NetworkConstruction.RData")


QTL <- read.table("eQTLinAllGenes.txt",header=T,row.names=1)
cauc.specific.egene.col <- labels2colors(QTL$Caucasian.specific.gene)
cns.specific.egene.col <- labels2colors(QTL$Chinese.specific.gene)
share.egene.col <- labels2colors(QTL$Shared)

pdf("Signed module colors for All gene dendrogram.pdf",height=7,width=11)
plotDendroAndColors(consTree, colors = cbind(moduleColors,cauc.signed.moduleColors, cns.signed.moduleColors,cauc.specific.egene.col,cns.specific.egene.col,share.egene.col), groupLabels = c("Consensus modules","Caucasian modules", "Chinese modules","Caucasian.specific.eGenes","Chinese.specific.eGenes","Overlap eGenes"), cex.colorLabels=0.6, cex.dendroLabels=0.2)
dev.off()

nSets = 2
setLabels = c("Caucasian", "Chinese")
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);
MET = consensusOrderMEs(consMEsC);
kME = consensusKME(multiExpr, universalColors = moduleColors)
pdf(file = "EigengeneNetworks.pdf", width= 8, height = 10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off();
