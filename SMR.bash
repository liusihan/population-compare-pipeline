for i in {1..5}; do ~/softwares/smr_Linux --eqtl-summary v"$i".eQTL --fastqtl-nominal-format --make-besd --out v"$i"; done

cat v1.esi | awk 'NR==FNR{a[$2]=0;b[substr($2,length($2))]=0}NR>FNR&&($2 in a){if($3 in b){split($2,c,"_");print $1"\t"$2"\t"0"\t"c[2]"\t"c[4]"\t"c[3]"\t"$5}else{split($2,c,"_");print $1"\t"$2"\t"0"\t"c[2]"\t"c[4]"\t"c[3]"\t"(1-$5)}}' - snp.frq >v1.esi.2
cat v1.epi | awk 'NR==FNR{a[$2]=0}NR>FNR&&($5 in a){print substr($1,4)"\t"$5"\t"0"\t"$2"\t"$8"\t"$4}' - ~/softwares/QTLtools/annotation.gene.gencodeV19.txt >v1.epi.2
cp v1.epi.2 v1.epi
cp v1.esi.2 v1.esi
~/softwares/smr_Linux --beqtl-summary v1 --update-epi v1.epi
~/softwares/smr_Linux --beqtl-summary v1 --update-esi v1.esi