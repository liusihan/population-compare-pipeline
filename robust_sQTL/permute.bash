#!/bin/sh
# MAKE SURE FUSION.compute_weights.R IS IN YOUR PATH
# FILL IN THESE PATHS
if [ $# -gt 5 ]; then
    echo "ERROR: YOU MUST INPUTE 5 PARAMETER like: bash permute.bash sample_size index PEER_factor FDR output_dir"
else
    ../permutation/index.R $1 $2
	mkdir $5
    cat qqnorm_all | awk 'NR==FNR&&NR>=1{a[$1]=0}NR>FNR{b=$1;for(i=2;i<=NF;i++){if((i-1) in a){b=b"\t"$i}};print b}' $2 - >$5/qqnorm_all
    ../permutation/RNA_processing.R $5/qqnorm_all $3 $5
	sed -n '1p' $5/qqnorm_all > $5/samplename
	cat $5/samplename $5/expr"$3".txt > $5/Covariant.txt
    awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[3]}' $5/qqnorm_all > $5/title.txt
	sed '1a #chr\tstart\tend' $5/title.txt |sed '1d' > $5/title
	paste -d "\t" $5/title $5/qqnorm_all > $5/qqnorm_all.bed
    for i in {1..22}; do cat $5/qqnorm_all.bed | awk 'NR==1{print $0}NR>=2&&$1=="'$i'"{print $0}' | sort -k1,1n -k2,2n | awk '$3>$2{$4=$4" . +"; print $0 }$2>$3{$4=$4" . -"; print $0}' | tr " " "\t" | bgzip >$5/phenotype.chr"$i".bed.gz;tabix -p bed $5/phenotype.chr$i.bed.gz;done
    for i in {1..22}; do     awk '          
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
    ' <(zcat $5/phenotype.chr1.bed.gz | awk 'NR==1{for(i=7;i<=NF;i++){b=b"""\t"""$i};print b}') <(zcat /zs32/data-analysis/liucy_group/shareData/BrainGVEX/eQTL/QTLtools/Genotype/genotypes.all.chr"$i".vcf.gz) | bgzip>$5/genotypes.all.chr$i.vcf.gz; tabix -p vcf $5/genotypes.all.chr$i.vcf.gz; done
    parallel -j 10 QTLtools cis --vcf $5/genotypes.all.chr{}.vcf.gz --bed $5/phenotype.chr{}.bed.gz --region {} --cov $5/Covariant.txt --out $5/sqtl.permute.chr{} --permute 1000 ::: {1..22}
    parallel -j 10 QTLtools cis --vcf $5/genotypes.all.chr{}.vcf.gz --bed $5/phenotype.chr{}.bed.gz --region {} --cov $5/Covariant.txt --out $5/sqtl.nopermute.chr{} --nominal 1 ::: {1..22}
    cat $5/sqtl.permute.chr* | gzip -c > $5/permutations_all.txt.gz
    /zs32/home/frwang/software/FastQTL/chinese_output_5/runFDR_cis.R $5/permutations_all.txt.gz $4 $5/permutations_all
    parallel -j 10 QTLtools cis --vcf $5/genotypes.all.chr{}.vcf.gz --bed $5/phenotype.chr{}.bed.gz --region {} --cov $5/Covariant.txt --out $5/sqtl.conditional.chr{} --mapping $5/permutations_all.thresholds.txt ::: {1..22} 
    cd $5
    ../FDR.R
    cat sqtl.nopermute.chr*.qvalue | awk '$NF<0.05' >sqtl.nopermute.all.signicifant
    cat sqtl.conditional.chr* | awk '$19==1' >sqtl.conditional.significant
    cd ..
fi
