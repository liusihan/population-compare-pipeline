smr_Linux --eqtl-summary Chinese.eQTL --qtltools-nominal-format --make-besd --out Chinese

cat Chinese.esi | awk 'NR==FNR{a[$2]=0;b[substr($2,length($2))]=0}NR>FNR&&($2 in a){if($3 in b){split($2,c,"_");print $1"\t"$2"\t"0"\t"c[2]"\t"c[4]"\t"c[3]"\t"$5}else{split($2,c,"_");print $1"\t"$2"\t"0"\t"c[2]"\t"c[4]"\t"c[3]"\t"(1-$5)}}' - snp.frq >Chinese.esi.2
cat Chinese.epi | awk 'NR==FNR{a[$2]=0}NR>FNR&&($5 in a){print substr($1,4)"\t"$5"\t"0"\t"$2"\t"$8"\t"$4}' - annotation.gene.gencodeV19.txt >Chinese.epi.2
cp Chinese.epi.2 Chinese.epi
cp Chinese.esi.2 Chinese.esi
smr_Linux --beqtl-summary Chinese --update-epi Chinese.epi
smr_Linux --beqtl-summary Chinese --update-esi Chinese.esi
