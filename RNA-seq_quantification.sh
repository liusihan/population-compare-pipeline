##merge fq.gz files and quality control
for i in RNA-seq/*
cd $i
zcat CL*_1.fq.gz | gzip > merge_1.fq.gz
zcat CL*_2.fq.gz | gzip > merge_2.fq.gz
/opt/tools/FastQC/fastqc -t 8 ./merge_*.fq.gz
cd ..
done


##STAR & RNA-seQC & picardtools:collect sequncing statistic metrics

for i in sample1 sample2 sample3; do /opt/tools/seq-analysis/STAR-2.5.2b/bin/Linux_x86_64/STAR --runThreadN 10 --genomeDir rsem/ref_hg19 --readFilesCommand zcat --readFilesIn $i/merge_1.fq.gz $i/merge_2.fq.gz --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA-seQC/$i --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20 && rm RNA-seQC/*.sam RNA-seQC/*.junction && rm RNA-seQC/*.toTranscriptome.out.bam && java -jar picard-tools-1.119/AddOrReplaceReadGroups.jar I=RNA-seQC/"$i"Aligned.sortedByCoord.out.bam O=RNA-seQC/"$i"Aligned.sortedByCoord.out_Pit.bam PL=illumine ID="$i"_sequence LB="$i"_sequence SM=Si PU=HWI-ST303 && java -jar picard-tools-1.119/ReorderSam.jar I=RNA-seQC/"$i"Aligned.sortedByCoord.out_Pit.bam O=RNA-seQC/"$i"Aligned.sortedByCoord.out_Pit_reorder.bam R=Reflib/ucsc.hg19.fasta CREATE_INDEX=TRUE && java -jar picard-tools-1.119/MarkDuplicates.jar INPUT=RNA-seQC/"$i"Aligned.sortedByCoord.out_Pit_reorder.bam OUTPUT=RNA-seQC/"$i"Aligned.sortedByCoord.out_Pit_re_marked.bam REMOVE_DUPLICATES=false METRICS_FILE=RNA-seQC/"$i".txt ASSUME_SORTED=true && samtools index RNA-seQC/"$i"Aligned.sortedByCoord.out_Pit_re_marked.bam; done


##RSEM: quanlity gene and isoform expression
rsem-calculate-expression --paired-end --star --star-path STAR-STAR_2.4.2a/STAR/bin/Linux_x86_64 --star-gzipped-read-file -p 8 merge_1.fq.gz merge_2.fq.gz ref_hg19/human "$i"_paired_end_quals

##combine sample results:awk
cat *.genes.results | awk -F '.' '{print $1}' | while read sample
>do 
>cat "$sample".genes.results | awk 'NR==1{print "'$sample'"}NR>=2{print $5}' >"$sample".count
>cat "$sample".genes.results | awk 'NR==1{print "'$sample'"}NR>=2{print $6}' >"$sample".tpm
>cat "$sample".genes.results | awk 'NR==1{print "'$sample'"}NR>=2{print $7}' >"$sample".fpkm
>done

cat sample1.genes.results | awk '{print $1}' >gene_ID
paste gene_ID *.count >raw_count
paste gene_ID *.tpm >raw_tpm
paste gene_ID *.fpkm >raw_fpkm
