##make directories 
mkdir corhot1
mkdir corhot2


#Make gene.tab files
cat Gene.ID | parallel -j 10 "cat EUR.nopermute.all | grep {} | awk 'BEGIN{print \"SNP A1 A2 P Effect N\"}\$1==\"{}\"{split(\$8,a,\"_\");print \$8,a[3],a[4],\$12,\$13,397}' > cohort2/{}.tab"

cat Gene.ID | parallel -j 10 "cat EAS.nopermute.all | grep {} | awk 'BEGIN{print \"SNP A1 A2 P Effect N\"}\$1==\"{}\"{split(\$8,a,\"_\");print \$8,a[3],a[4],\$12,\$13,145}' > cohort1/{}.tab"

#perform meta-analysis with METAL
for i in {1..22}
do
cat Gene_position.txt | awk '$2=='$i'{print $1}'  | while read gene; do ln -s -f cohort1/"$gene".tab cohort1.stat; ln -s -f cohort2/"$gene".tab cohort2.stat; ../generic-metal/executables/metal metal.txt >metal.log; cat METAANALYSIS1.TBL | awk '{print "'$gene'",$0}' >>meta_eqtl.chr"$i"; done
done

#FDR
library(qvalue)
for(i in c(1:22)){
    d = read.table(paste("eqtl.permute.chr",i,sep=""), head=F, sep=" ", stringsAsFactors=F)
    d$qvalue = p.adjust(d$V11,method = "fdr")
    write.table(d[which(!is.na(d$qvalue)),], paste("eqtl.permute.chr",i,".qvalue",sep=""), quote=F, row.names=F, col.names=T, sep="\t")
}