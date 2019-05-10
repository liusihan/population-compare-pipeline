#!/bin/sh
# MAKE SURE index.r RNA_processing.r and FDR.r ARE IN YOUR PATH
# FILL IN THESE PATHS
if [ $# -ne 5 ]; then
    echo "ERROR: YOU MUST INPUTE 5 PARAMETER like: bash permute.bash sample_size index.txt PEER_factor FDR output_dir raw_counts.txt"
else
    Rscript index.r $1 $2
    mkdir $5
    cat $6 | awk 'NR==FNR&&NR>=1{a[$1]=0}NR>FNR{b=$1;for(i=2;i<=NF;i++){if((i-1) in a){b=b"\t"$i}};print b}' $2 - >$5/expression_raw_count
    Rscript RNA_processing.r $5/expression_raw_count $3 $5
    cat $5/expr"$3".xls | awk 'BEGIN{OFS="\t"}NR==1{$1="ID""\t"$1;print $0}NR>=2{print $0}' >$5/Covariant.txt
    awk 'BEGIN{FS="\t";OFS="\t"}NR==FNR{a[$1]=1}NR>FNR{if($4 in a){start=$3;$3=$2;$2=start}; print $0}' <(zcat gencode.v19.annotation.bed.gz |awk 'BEGIN{FS="[\t.]"}$5=="gene"{a[$7]=$6}END{for(i in a){if(a[i]=="-") print i}}') <(awk 'BEGIN{FS="\t";OFS="\t"}NR==FNR&&$5=="gene"{split($8,s,".");a[s[1]]=$1"\t"$2"\t"$3}NR>FNR&&FNR==1{print "#chr\tstart\tend\tgene"$0}NR>FNR&&FNR>1&&($1 in a){print a[$1],$0}' <(zcat gencode.v19.annotation.bed.gz) $5/log2cpm.fgene.fsample.qn) |sed -e's/^chr//'>$5/log2cpm.fgene.fsample.qn.bed
    for i in {1..22} X; do cat $5/log2cpm.fgene.fsample.qn.bed | awk 'NR==1{print $0}NR>=2&&$1=="'$i'"{print $0}' | sort -k1,1n -k2,2n | awk '$3>$2{$4=$4" . +"; print $0 }$2>$3{$4=$4" . -"; print $0}' | tr " " "\t" | bgzip >$5/phenotype.chr"$i".bed.gz;tabix -p bed $5/phenotype.chr$i.bed.gz;done
    for i in {1..22} X; do     awk '          
        BEGIN{FS="\t";OFS="\t"}
        NR==FNR{for(i=2;i<=NF;i++) order[$i]=i+8}
        NR>FNR{
            if($1~/^##/){
                print $0
            }else if($1=="#CHROM"){
                for(i=10;i<=NF;i++){
                    if($i in order){
                        ind[order[$i]]=i
                    }
                }
                printf($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9)
                for(i=10;i<=length(ind)+9;i++){
                    n=ind[i]
                    printf("\t"$n)
                }
                printf("\n")
            }else{
                if($5~/,/) next
                $3=$1"_"$2"_"$4"_"$5;
                printf($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9)
                for(i=10;i<=length(ind)+9;i++){
                    n=ind[i]
                    printf("\t"$n)
                }
                printf("\n")
            }
        }
    ' <(zcat $5/phenotype.chr1.bed.gz | awk 'NR==1{for(i=7;i<=NF;i++){b=b"""\t"""$i};print b}') <(zcat genotypes.all.chr"$i".vcf.gz) | bgzip>$5/genotypes.all.chr$i.vcf.gz; tabix -p vcf $5/genotypes.all.chr$i.vcf.gz; done
    parallel -j 10 QTLtools_1.0_CentOS6.8_x86_64 cis --vcf $5/genotypes.all.chr{}.vcf.gz --bed $5/phenotype.chr{}.bed.gz --region {} --cov $5/Covariant.txt --out $5/eqtl.permute.chr{} --permute 1000 ::: {1..22} X
    parallel -j 10 QTLtools_1.0_CentOS6.8_x86_64 cis --vcf $5/genotypes.all.chr{}.vcf.gz --bed $5/phenotype.chr{}.bed.gz --region {} --cov $5/Covariant.txt --out $5/eqtl.nopermute.chr{} --nominal 1 ::: {1..22} X
    cat $5/eqtl.permute.chr* | gzip -c > $5/permutations_all.txt.gz
    Rscript QTLtools/script/runFDR_cis.R $5/permutations_all.txt.gz $4 $5/permutations_all
    parallel -j 10 QTLtools_1.0_CentOS6.8_x86_64 cis --vcf $5/genotypes.all.chr{}.vcf.gz --bed $5/phenotype.chr{}.bed.gz --region {} --cov $5/Covariant.txt --out $5/eqtl.conditional.chr{} --mapping $5/permutations_all.thresholds.txt ::: {1..22} X
    cd $5
    Rscript ../FDR.r
    cat eqtl.nopermute.chr*.qvalue | awk '$NF<0.05' >eqtl.nopermute.all.signicifant
    cat eqtl.conditional.chr* | awk '$19==1' >eqtl.conditional.significant
    cd ..
fi
