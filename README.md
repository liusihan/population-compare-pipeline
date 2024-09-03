# Cross-ancestry analysis of brain QTLs enhances interpretation of schizophrenia GWAS

This repository contains analysis and ploting code in our paper about constructing brain expression regulation architecture in African American, European, and East Asian populations:<br>
  * DNA-seq quality control and imputation<br>
  * RNA-seq alignment, quantification, and quality control<br>
  * eQTL and sQTL mapping<br>
  * Intergrating with GWAS summary result<br>
  * Preservation test and robust WGCNA<br>

## Introduction
Research on brain expression quantitative trait loci (eQTLs) has illuminated the genetic underpinnings of schizophrenia (SCZ). Yet, most of these studies have been centered on European populations, leading to a constrained understanding of population diversities and disease risks. To address this gap, we examined genotype and RNA-seq data from African Americans (AA, n=158), Europeans (EUR, n=408), and East Asians (EAS, n=217). The two key questions we sought to answer are: 1) what drives the brain eQTL differences across populations? 2) what do we gain by studying brain eQTLs in diverse populations?

## Install

```Linux
git clone https://github.com/liusihan/EAS-and-EUR-brain-regulatory-pattern-discovory-and-comparson-pipeline.git
```

## Dependencies
  * R-3.3.3
  * [QTLtools](https://qtltools.github.io/qtltools/)
  * [LDSC](https://github.com/bulik/ldsc)
  * [SMR](https://cnsgenomics.com/software/smr/)
  * [PrediXcan](https://github.com/hakyim/PrediXcan)


## Usage

### QTL replication rate
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

## Data Availability
The summary statistics of eQTLs are available on the CNBrainQTL (Chinese Brain xQTL) database at http://121.41.22.25:8080/qtl_v1/#/ or can be downloaded from the `eQTL_results` folder.

