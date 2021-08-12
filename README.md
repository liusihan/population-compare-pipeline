# Analysis pipelines for eQTL/sQTL discovery and GWAS signals interpratation

This repository contains analysis and ploting code in our paper about constructing brain expression architecture in East Asian population:<br>
  * DNA-seq quality control and imputation<br>
  * RNA-seq alignment, quantification, and quality control<br>
  * eQTL and sQTL mapping<br>
  * Intergrating with GWAS summary result<br>
  * Preservation test and robust WGCNA<br>

## Introduction

## Install

```Linux
git clone https://github.com/liusihan/EAS-and-EUR-brain-regulatory-pattern-discovory-and-comparson-pipeline
```

## Dependencies
  * R-3.3.3
  * [QTLtools](https://qtltools.github.io/qtltools/)
  * [LDSC](https://github.com/bulik/ldsc)
  * [SMR](https://cnsgenomics.com/software/smr/)
  * [PrediXcan](https://github.com/hakyim/PrediXcan)


## Usage

### QTL replicatin rate
```R
Rscript pi1.r P.txt
```

`NOTE`: P.txt is a text file with nominal P value per line without header. 


### Downsampling analysis
```Linux
bash robust_eQTL.sh sample_size index.txt PEER_factor FDR output_dir raw_counts.txt
```
* sample_size: the sample size you want to perform eQTL analysis
* index.txt: output file for sample index
* PEER_factor: number of peer hidden factor to collect
* FDR: storey FDR qvalue for adjust P value
* output_dir: prefix of output directory
* raw_counts.txt: raw counts data, each column stand for a sample, each raw stand for a gene

`NOTE`: the output_file include three type eQTL results under different mode in QTLtools(nominal/permutation/conditional)
