#Add an annotation column(eSNPs or sSNPs) in baseline model
for i in {1..22}; do zcat ~/database/LDSR/EUR_baseline_2.2/baselineLD."$i".annot.gz | awk 'BEGIN{OFS="\t"}NR==FNR{a[$9""$10]=0}NR>FNR&&FNR==1{print $0,"robust_eQTL"}NR>FNR&&FNR>=2{if($1""$2 in a){print $0,1}else{print $0,0}}' nominal_eQTL.txt - | gzip >~/database/LDSR/my_baseline_EUR/baselineLD."$i".annot.gz; done

#Calculated LD score
for i in {1..22}; do python ~/ldsc/ldsc.py --l2 --bfile ~/database/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC."$i" --ld-wind-cm 1 --annot ~/database/LDSR/my_baseline_EUR/baselineLD."$i".annot.gz --print-snps ~/database/LDSR/hapmap3_snps/hm."$i".snp --out ~/database/LDSR/my_baseline_EUR/baselineLD."$i"; done

#Make sumstats file
python ~/ldsc/munge_sumstats.py --sumstats PGC3.SCZ.summary --N 65967 --out PGC3.SCZ --merge-alleles ~/database/LDSR/w_hm3.snplist --a1-inc

#Partitioned heritability
python ~/ldsc/ldsc.py --h2 PGC3.SCZ.sumstats.gz --ref-ld-chr ~/database/LDSR/my_baseline_EUR/baselineLD. --w-ld-chr /zs32/home/shliu/database/LDSR/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr ~/database/LDSR/1000G_Phase3_frq/1000G.EUR.QC. --out PGC_SCZ
