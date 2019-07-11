library(WGCNA)

load("BrainGVEx_human_brain_AdjustCombatandCovariants.rdata")
load("../tpmoverlapgenes.rdata")

datExpr <- AA_datExpr[genes,]

# Choose a set of soft-thresholding powers
powers = c(c(1:30), seq(from = 20, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(t(datExpr), powerVector = powers, networkType = "unsigned",verbose = 5)
# Plot the results:
pdf("./unsigned_softThresh.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

net = blockwiseModules(t(datExpr), power = 6,mergeCutHeight=0.25,nThreads=10,
                                maxBlockSize=50000,networkCalibration  = 'single quantile',
                                calibrationQuantile = 0.8,numericLabels = TRUE,
                                corType="bicor", useMean = T,
                                minModuleSize=20,pamStage=FALSE,networkType="signed",saveTOMs = TRUE,saveTOMFileBase = "./CAUCBrainSignedTOM",verbose = 3)



load("CNSBrainUnSignedTOM-block.1.RData")
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
colors = vector(mode="list")
labels = vector(mode="list")

for (pam in c(TRUE,FALSE)) {
  for (minModSize in c(20,30,40,50)) {
    for (dthresh in c(0.1, 0.2)) {
      for(ds in c(1,2,3,4)) { 
        print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-TOM))
        merged = mergeCloseModules(exprData= t(datExpr), colors = tree$labels, cutHeight=dthresh)
        colors = cbind(colors, labels2colors(merged$colors))
        
        labels = c(labels, paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
      }
    }
  }
}
pdf("Fig2-UnSigned WGCNA different params.pdf", wi = 15)
plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE, dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.5)
dev.off()


colors = vector(mode="list")
labels = vector(mode="list")

minModSize = 100
pamStage = F
ds = 4
dthresh = 0.3

tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pamStage, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-TOM)); 
merged = mergeCloseModules(exprData= t(datExpr), colors = tree$labels, cutHeight=dthresh)
colors = cbind(colors, labels2colors(merged$colors))
colors = labels2colors(merged$colors)
table(colors)
length(table(colors))

pdf("UnSigned Dendrogram and the Module Colors.pdf")
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(geneTree, colors,
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = merged$colors
moduleColors = labels2colors(merged$colors)
MEs = merged$newMEs;
geneTree = geneTree;
save(minModSize,pamStage,ds ,dthresh, MEs, moduleLabels, moduleColors, geneTree,merged,file = "CNS-UnSigned-networkConstruction-auto.RData")
