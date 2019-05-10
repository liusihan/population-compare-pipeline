for k in 5000 10000 20000 500000 1000000; do
  cat /zs32/data-analysis/liucy_group/shareData/Chinese_brain/eQTL_result/QTLtools/result/conditional/Chinese_sig_specific.txt | awk '{if($4>$3){print $2,$3-'$k',$4+'$k'}else{print $2,$4-'$k',$3+'$k'}}' >Chinese_specific_eGene_"$k".bed;
  for i in {1..22}; do cat ~/database/LDSR/1000G_Phase3_EAS_baseline_v1.2_ldscores/"$i".txt | awk 'NR==FNR&&$1=='$i'{a[NR]=$2;b[NR]=$3}NR>FNR&&FNR==1{print "Chinese"}NR>FNR&&FNR>=2{flag=0;for(i in a){if(a[i]<$1&&b[i]>$1){flag=1}};print flag}' Chinese_specific_eGene_"$k".bed - >chinese_annotate_"$i".txt; done;
  for i in {1..22}; do zcat ~/database/LDSR/1000G_Phase3_EAS_baseline_v1.2_ldscores/baseline."$i".annot.gz | awk 'BEGIN{OFS="\t"}NR==FNR{a[NR]=$0}NR>FNR{print a[FNR],$1}' - chinese_annotate_"$i".txt | gzip >baseline"$k"."$i".annot.gz; done;
  for i in {1..22}; do python ~/ldsc/ldsc.py --l2 --bfile ~/database/LDSR/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC."$i" --ld-wind-cm 1 --annot baseline"$k"."$i".annot.gz --print-snps ~/database/LDSR/hapmap3_snps/hm."$i".snp --out baseline"$k"."$i"; done;  
  python ~/ldsc/ldsc.py --h2 ../PGC.eas.SCZ.sumstats.gz --ref-ld-chr baseline"$k". --w-ld-chr /zs32/home/shliu/database/LDSR/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC. --overlap-annot --frqfile-chr ~/database/LDSR/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC. --out EAS.SCZ_"$k"; done

