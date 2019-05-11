Chinese-brain-project eQTL pipeline
====

This repository contains analysis pipelines for:<br>
  * RNA-seq alignment, quantification, and quality control<br>
  * eQTL mapping and functional annotation<br>
  * Intergrate with GWAS summary result<br>
  * Preservation test and robust WGCNA<br>


### Install
```Linux
git clone https://github.com/liusihan/Chinese-brain-project
```

## Usage

### QTl replicatin rate
```R
Rscript pi1.r P.txt
```
NOTE: P.txt is a text file with nominal P value per line without head


### robust eQTL
```Linux
bash robust_eQTL.sh sample_size index.txt PEER_factor FDR output_dir raw_counts.txt
```
* sample_size: the sample size you want to perform eQTL analysis
* index.txt: output file for sample index
* PEER_factor: number of peer hidden factor to collect
* FDR: storey FDR qvalue for adjust P value
* output_dir: prefix of output directory
* raw_counts.txt: raw counts data, each column stand for a sample, each raw stand for a gene
  
