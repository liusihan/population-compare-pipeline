#3g_Annotation_Pathway.R

rm(list=ls()); options(stringsAsFactors = F)
library(pSI); library(gProfileR); library(gplots); library(biomaRt); library(WGCNA); library(stringr); library(GENIE3)

load("Chinese-Signed-networkConstruction-auto.RData")
load("Chinese_human_brain_Adjust_Expr&Meta.rdata")
load("BrainGVEx_human_brain_AdjustCombatandCovariants.rdata")
load("../overlapgenes.rdata")
datExpr <- datExpr[genes,]

eigmat = MEs; colnames(eigmat) = gsub("ME","",colnames(eigmat))
kME = signedKME(t(datExpr), MEs); colnames(kME) = gsub("kME", "", colnames(kME))

ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
bm = getBM(attributes = c("ensembl_gene_id", "external_gene_id"), filters = "ensembl_gene_id", values = a$gene, mart=ensembl)
datProbes = data.frame()
datProbes = data.frame(ensg=rownames(datExpr))
datProbes = data.frame(ensg=rownames(datExpr), symbol=bm$external_gene_id[match(rownames(datExpr), bm$ensembl_gene_id)])
rownames(datProbes) <- rownames(datExpr)

## Calculate GO enrichment of Race-associated modules
resultsGO = data.frame()
for (m in c("pink" )) {
  me_name = paste("ME", m, sep="")
  i = which(m == unique(moduleColors))  
  
  ## GO Enrichment
  query = rownames(datExpr)[moduleColors==m]
  go = gprofiler(query, organism="hsapiens", ordered_query = F, significant = T, exclude_iea = F, region_query = F,max_p_value = 1,correction_method = "fdr", custom_bg = rownames(datExpr),
                 max_set_size = 2000,hier_filtering = "moderate", domain_size = "annotated", 
                 numeric_ns = "", include_graph = F,src_filter = c("GO", "KEGG"))
  go = go[order(go$p.value),]
  
  num_to_plot=10
  par(oma=c(0,10,0,0))
  bp = barplot(-log10(as.numeric(go$p.value[num_to_plot:1])), main=paste("\n\n\n\n\n",m, "- GO, KEGG"), horiz=T, yaxt='n', col=m, xlab='-log10(FDR)\n',cex.main=0.7, cex.axis = .7)
  axis(2,at=bp,labels=str_wrap(go$term.name[num_to_plot:1]),tick=FALSE,las=2,cex.axis=.8);
  abline(v=-log10(0.01), col="red", lwd=2,lty=2)
  
  resultsGO = rbind(resultsGO, cbind(m, go))
}
write.csv(file="../NetworkComparision/ModuleTrait/Caucasian-Signed-network_gProfiler_resultsGO022619.csv", resultsGO)

resultsGOsub2 = data.frame()
for(c in sort(c("pink" ))) {
  idx = which(resultsGO$m==c) 
  resultsGOsub2 = rbind(resultsGO[idx[5:1],],resultsGOsub2)  
}


pdf("../NetworkComparision/ModuleTrait/Chinese-Signed-specificity-Module-GO-pink.pdf",width = 4,height=2)
par(oma=c(0,6,0,0), mar=c(5,4,2,2))
bp = barplot(-log10(resultsGOsub2$p.value), main="Chinese-Signed-Pathway-Enrichment", horiz=T, yaxt='n', col=resultsGOsub2$m, xlab='-log10(p)',cex.main=0.5, cex.axis =0.5, border=NA)
axis(2,at=bp,labels=str_wrap(resultsGOsub2$term.name,width=45),tick=FALSE,las=2,cex.axis=0.5);
abline(v=-log10(0.01), col="black", lwd=2,lty=2)
dev.off()


# Calculate Cell-Type specificity of modules
#Zhang using pSI
zhang.datExpr = read.csv("../Annotation/datExpr.zhangHuman.avgForPSI.csv",row.names=1)[,-1]
set.seed(100)
pSI.output = specificity.index(pSI.in=zhang.datExpr,bts=100,p_max=.1, e_min=0.3); 
pSI.count(pSI.output)

cell.p.zhang = matrix(NA, length(unique(moduleColors)),5);  rownames(cell.p.zhang) = unique(moduleColors)
colnames(cell.p.zhang) = colnames(pSI.output)


for(mod in rownames(cell.p.zhang)) {
    f = fisher.iteration(pSI.output, rownames(datExpr)[moduleColors==mod],p.adjust = F)
    cell.p.zhang[mod,] = f$`0.05 - nominal`
}

cell.p.zhang.fdr = p.adjust(cell.p.zhang,"fdr")
dim(cell.p.zhang.fdr) = dim(cell.p.zhang); dimnames(cell.p.zhang.fdr) = dimnames(cell.p.zhang);
to_plot = cell.p.zhang.fdr[c("black" ,"magenta", "tan" ),]

colnames(eigmat) <- labels2colors(names(eigmat))
dendro.col = as.dendrogram(hclust(as.dist(1-bicor(zhang.datExpr)), method="average"))
denro.row= as.dendrogram(hclust(as.dist(1-bicor(eigmat[,c("black" ,"magenta", "tan" )])),method="average"))

pdf("../NetworkComparision/ModuleTrait/Caucasian-Signed-specificity-Module-CellType.pdf",width=6,height=5)
heatmap.2(-log10(to_plot),col=blueWhiteRed(1000,1)[500:1000],
          scale="none",trace="none",cexRow = 0.8,cexCol = .8, density.info = "none",
          colsep=0:7,rowsep=0:8,sepcolor="grey",sepwidth=c(0.02,0.02),
          srtCol=45,offsetRow=0,offsetCol=-0.5,
          Rowv=denro.row, Colv=dendro.col,
          key=T,key.xlab="-log10(P)", cellnote=signif(to_plot,1), notecex=.8, notecol="black",main="Caucasian-Signed-Enrichment")
dev.off()


out = cell.p.zhang.fdr
colnames(out) = paste("Enrichment.", colnames(out), ".FDR",sep="")
write.csv(out,file="../NetworkComparision/ModuleTrait/Caucasian-Signed-CellType.csv")

## Calculate TFBS enrichment of Disease-associated modules
resultsTFBS = data.frame()
for (m in c("black" ,"magenta", "tan" )) {
  me_name = paste("ME", m, sep="")
  i = which(m == unique(moduleColors))  
  
  ## GO Enrichment
  query = rownames(datExpr)[moduleColors==m]
  go = gprofiler(query, organism="hsapiens", ordered_query = F, significant = T, exclude_iea = F, region_query = F,max_p_value = 1,correction_method = "fdr", custom_bg = rownames(datExpr),
                 hier_filtering = "strong", domain_size = "annotated", numeric_ns = "", include_graph = F,src_filter = c("TF"))
  go = go[order(go$p.value),]
  
  resultsTFBS = rbind(resultsTFBS, cbind(m, go[1:min(nrow(go),20),]))
}
write.csv(resultsTFBS,file="../NetworkComparision/ModuleTrait/Caucasian-Signed-TableS2-TFBS.csv")


colnames(kME) <- labels2colors(names(kME))
## Identify hub genes transcription factors
hubGenes = data.frame(Module=NA, Gene=NA, Symbol=NA,Rank=NA)
for (m in c("black" ,"magenta", "tan" )) {
  mod = rownames(datExpr)[moduleColors==m]
  mod = mod[order(kME[mod,m],decreasing = T)[1:20]]
  sym = datProbes[mod,2]
  hubGenes=rbind(hubGenes, data.frame(Module=m, Gene=mod, Symbol = sym, Rank=1:20))
}
hubGenes = hubGenes[-1,]


ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
a=listAttributes(ensembl);f=listFilters(ensembl)
bm = getBM(attributes = c("ensembl_gene_id", "go_id", "go_linkage_type","goslim_goa_accession","goslim_goa_description"), filters = "ensembl_gene_id", values = hubGenes$Gene, mart=ensembl)
bm = bm[grep("transcription factor",bm$goslim_goa_description),]
hubGenes$TF = bm$goslim_goa_description[match(hubGenes$Gene, bm$ensembl_gene_id)]

write.csv(hubGenes,file="../NetworkComparision/ModuleTrait/Caucasian-Signed-HubGeneTFs.csv", row.names = F)

# incorporates TF information to construct a regulatory network
moduleLink = data.frame()
for (m in c("black" ,"magenta", "tan" )) {
  me_name = paste("ME", m, sep="")
  i = which(m == unique(moduleColors))  
  
  ## GENIE3 Enrichment
  query = rownames(datExpr)[moduleColors==m]
  genieExpr = datExpr[query,]
  set.seed(123)
  weightMat <- GENIE3(genieExpr, regulators = intersect( hubGenes$Gene[which(hubGenes$Module==m)], hubGenes$Gene[which(hubGenes$TF!="NA")]), verbose=TRUE)
  linkList <- getLinkList(weightMat, threshold=0.3)
  moduleLink <- rbind(moduleLink,data.frame(linkList,module=m))
}
moduleLink <- data.frame(moduleLink,regulatory.symbol=datProbes$symbol[match(moduleLink$regulatoryGene,datProbes$ensg)],target.symbol=datProbes$symbol[match(moduleLink$targetGene,datProbes$ensg)])
write.csv(moduleLink,file="../NetworkComparision/ModuleTrait/Caucasian-Signed-DirectNetwork.csv", row.names = T)

#disease genes enrichment
eRNAnetwork <- read.csv("../Annotation/DiseaseGeneList.csv",header=T)
source("../Annotation/fisher_overlap.R")
table.p = matrix(NA, nrow=1, ncol=dim(eRNAnetwork)[2])
rownames(table.p) = "pink"; 
colnames(table.p) = colnames(eRNAnetwork)
table.or = table.p
hgnc = rownames(datExpr)
m="pink"
for(e in (1:ncol(table.p))) {
	x = length(intersect(hgnc[moduleColors=="pink"],eRNAnetwork[,e]))
	d = length(hgnc)-length(eRNAnetwork[,e])
	n = length(eRNAnetwork[,e][which(eRNAnetwork[,e]!="")])
	k = length(hgnc[moduleColors==m])
    f = dhyper(x,n,d,k)
    table.or[m,e] = length(intersect(hgnc[moduleColors==m],eRNAnetwork[,e]))
    table.p[m,e] = as.numeric(f)
}


table.p.fdr = p.adjust(table.p,"fdr")
dim(table.p.fdr) = dim(table.p); dimnames(table.p.fdr) = dimnames(table.p)
table.p.fdr[table.or<1] = 1

to_plot = as.data.frame(table.p.fdr)

p=ggplot(melt(to_plot),aes(x=variable,y=-log10(value),fill="pink")) + 
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  ggtitle("Chinese Signed Disorder Genes Enrichment") + theme_classic() +
  geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
  theme(axis.text.x=element_text(angle=55, size=10, hjust=1),legend.position = "none")
p

ggsave(p, filename = "../NetworkComparision/ModuleTrait/Chinese-Signed-specificity-Module-DisorderGene-pink.pdf",width=12,height=5)
 
