##phenotype is prepared by leafcutter
for bamfile in $(ls ~/*.bam);do£» echo Converting $bamfile to $bamfile.junc; sh bam2junc.sh $bamfile $bamfile.junc; echo $bamfile.junc >> Chinese_juncfiles.txt; done
python /leafcutter/clustering/leafcutter_cluster.py -j Chinese_juncfiles.txt -m 50 -o ../output/Chinese -l 500000
python scripts/prepare_phenotype_table.py output/Chinese_perind.counts.gz -p 10

##sQTL mapping
#nominal
parallel -j 10 QTLtools cis --vcf genotypes.all.chr{}.vcf.gz --bed Chinese.chr{}.bed.gz --region {} --cov Covariant.txt --out eqtl.permute.chr{} --nominal 1 ::: {1..22}

#permutation
parallel -j 10 QTLtools cis --vcf genotypes.all.chr{}.vcf.gz --bed Chinese.chr{}.bed.gz --region {} --cov Covariant.txt --out eqtl.permute.chr{} --permute 1000 ::: {1..22}


