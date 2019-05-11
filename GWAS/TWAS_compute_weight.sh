#!/bin/sh
# MAKE SURE FUSION.compute_weights.R IS IN YOUR PATH
# FILL IN THESE PATHS
if [ $# -ne 3 ]; then
    echo "ERROR: YOU MUST INPUTE 3 PARAMETER like: bash TWAS_compute_weights.bash expression GENO covariant.txt"
else
    GCTA="gcta_nr_robust"
    PLINK="plink"
    PLINK2="plink2"
    GEMMA="gemma.linux"
# ALTERNATIVELY: ENSURE THAT plink, gcta, gemma CAN BE CALLED FROM PATH AND REMOVE --PATH_* FLAGS BELOW
# PATH TO DIRECTORY CONTAINING LDREF DATA (FROM FUSION WEBSITE or https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2)
    LDREF="TWAS/LDREF"
# THIS IS USED TO RESTRICT INPUT SNPS TO REFERENCE IDS ONLY

# PATH TO GEUVADIS GENE EXPRESSION MATRIX:
    PRE_GEXP=$1
# GEUVADIS DATA WAS DOWNLOADED FROM https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/

# PATH TO PREFIX FOR GEUVADIS GENOTYPES SPLIT BY CHROMOSOME
# SUBSAMPLE THESE TO THE LDREF SNPS FOR EFFICIENCY
    PRE_GENO=$2

# PATH TO OUTPUT DIRECTORY (population-specific subdirs will be made)


# ROWS IN THE MATRIX TO ANALYZE (FOR BATCHED RUNS)

# --- BEGIN SCRIPT:

# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:


# Loop through each gene expression phenotype in the batch
    cat $PRE_GEXP | awk '{print $0}' |  while read PARAM; do

# Get the gene positions +/- 500kb
    CHR=`echo $PARAM | awk '{ print $1 }'`
    P0=`echo $PARAM | awk '$2>0.5e6{ print $2 - 0.5e6 }$2<0.5e6{print 0}'`
    P1=`echo $PARAM | awk '{ print $3 + 0.5e6 }'`
    GNAME=`echo $PARAM | awk '{ print $4 }'`

    OUT="$PRE_GEXP.$GNAME"

    echo $GNAME $CHR $P0 $P1>>gene_pos.txt

# Pull out the current gene expression phenotype
    echo $PARAM | tr ' ' '\n' | tail -n+5 | paste $PRE_GEXP.ID - > $OUT.pheno

# Get the locus genotypes for all samples and set current gene expression as the phenotype
    $PLINK --bfile $PRE_GENO$CHR --allow-no-sex --pheno $OUT.pheno --make-bed --out $OUT --keep $OUT.pheno --chr $CHR --from-bp $P0 --to-bp $P1 --noweb --extract $LDREF/1000G.EUR.$CHR.bim

# Process all samples together (for reference purposes only since this is mult-ethnic data)
    Rscript FUSION.compute_weights2.R --bfile $OUT --PATH_plink $PLINK2 --tmp $OUT.tmp --out ./WEIGHTS/$OUT.final --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --covar $3 --models blup,lasso,top1,enet

# ALTERNATIVELY ADD COVARIATES HERE USING THE --covar FLAG
# MINIMAL COMMAND IS: `Rscript /zs32/home/shliu/softwares/fusion_twas-master/FUSION.compute_weights.R --bfile $OUT --tmp $OUT.$pop.tmp --out $FINAL_OUT`
# Remove all intermediate files
    rm $OUT.*
    done
fi


# GO TO THE NEXT GENE
