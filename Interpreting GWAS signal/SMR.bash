smr_Linux --eqtl-summary EAS.eQTL --qtltools-nominal-format --make-besd --out EAS

#snp.frq is a text file contains the frequency information which generated from plink
cat EAS.esi | awk 'NR==FNR{a[$2]=0;b[substr($2,length($2))]=0}NR>FNR&&($2 in a){if($3 in b){split($2,c,"_");print $1"\t"$2"\t"0"\t"c[2]"\t"c[4]"\t"c[3]"\t"$5}else{split($2,c,"_");print $1"\t"$2"\t"0"\t"c[2]"\t"c[4]"\t"c[3]"\t"(1-$5)}}' - snp.frq >EAS.esi.2
cat EAS.epi | awk 'NR==FNR{a[$5]=$8"""\t"""$4}NR>FNR&&($2 in a){print $1"\t"$2"\t"$3"\t"$4"\t"a[$2]}' ~/softwares/QTLtools/annotation.gene.gencodeV19.txt - |sort -n -k 1 -k 4 >EAS.epi.2
cp EAS.epi.2 EAS.epi
cp EAS.esi.2 EAS.esi
smr_Linux --beqtl-summary EAS --update-epi EAS.epi
smr_Linux --beqtl-summary EAS --update-esi EAS.esi
~/softwares/SMR/smr_Linux --bfile chrall --gwas-summary EAS_SCZ.ma --beqtl-summary EAS --out EAS_SCZ --thread-num 10
