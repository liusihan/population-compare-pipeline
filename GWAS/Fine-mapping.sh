cat /zs32/data-analysis/liucy_group/shareData/Chinese_brain/eQTL_result/Chinese_brain_significant.txt | awk '!a[$1]++' | grep -v X | awk '{print $1,$2}' | while read gene chr; do 
    cat /zs32/data-analysis/liucy_group/shareData/Chinese_brain/eQTL_result/Chinese_brain_significant.txt | awk '$1=="'$gene'"{print $8,$13}' >"$gene".txt; 
    Rscript Zscore.r "$gene".txt "$gene".Z.txt; 
    cat $gene.Z.txt | sed 's/"//g' >$gene.rs.Z.txt; cat $gene.rs.Z.txt | awk '{print $1}' >eSNP.ID; 
    plink --bfile chr"$chr" --extract eSNP.ID --make-bed --out $gene; 
    plink --bfile "$gene" --r2 --out "$gene"; 
    CAVIAR -l "$gene".ld -z "$gene".rs.Z.txt -o "$gene"; 
    cat "$gene"_set >>95credibleset.txt; 
    cat "$gene"_prop >>95credibleset.prop; 
    rm "$gene".*; 
    rm "$gene"_*; done
