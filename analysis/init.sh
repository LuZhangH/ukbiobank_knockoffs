#!/bin/bash

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22

echo "Initializing list of variants and individuals for BOLT-LMM"
echo "  First CHR         :" $CHR_MIN
echo "  Last CHR          :" $CHR_MAX

##################################
# Assemble covariate information #
##################################

source activate ukb
Rscript --vanilla init_covariates.R

#######################
# Initialize BOLT-LMM #
#######################

## Make list of SNPs to exlude ##

# Input
GENO_FILE="/scratch/groups/candes/ukbiobank/genotypes/ukb_gen_chr"

# Make directory for output
SNP_DIR="/scratch/groups/candes/ukbiobank_tmp/analysis/variants"
mkdir -p $SNP_DIR

# # Work on one chromosome at a time
# for CHR in $(seq $CHR_MIN $CHR_MAX); do
#   echo "Making list of excluded SNPs from chromosome $CHR ..."
#   # Remove knockoffs from the list of original variables used in this experiment
#   BIM_FILE="/scratch/groups/candes/ukbiobank_tmp/augmented_data/Radj50_K50/ukb_gen_chr"$CHR".bim"
#   SNP_FILE=$SNP_DIR"/ukb_gen_chr"$CHR".bim"
#   awk -F"[\t.]" '{print $2}' $BIM_FILE | uniq > $SNP_FILE
#   # List the names of the SNPs in the original data
#   SNP_TMP_FILE=$SNP_DIR"/ukb_gen_chr"$CHR"_tmp.txt"
#   awk -F"[\t.]" '{print $2}' $GENO_FILE$CHR".bim" > $SNP_TMP_FILE
#   # Find the names of the SNPs to exclude
#   SNP_EXCLUDE_FILE=$SNP_DIR"/ukb_gen_chr"$CHR"_exclude.txt"
#   comm -23 <(sort -i $SNP_TMP_FILE) <(sort -i $SNP_FILE) > $SNP_EXCLUDE_FILE
#   # Remove temporary files
#   rm $SNP_FILE $SNP_TMP_FILE
#   echo "List written on: "$SNP_EXCLUDE_FILE
# done

# Combine list of files for different chromosomes
SNP_EXCLUDE_FILE=$SNP_DIR"/ukb_gen_exclude.txt"
cat $SNP_DIR"/ukb_gen_chr"*"_exclude.txt" > $SNP_EXCLUDE_FILE
echo "Final list of exluded variants written on: "
echo $SNP_EXCLUDE_FILE

## Make list of individuals in the genotype files ##
mkdir -p "/scratch/groups/candes/ukbiobank_tmp/analysis/individuals"
FAM_FILE="/scratch/groups/candes/ukbiobank_tmp/analysis/individuals/ukb_individuals.fam"
awk '{print $1,$2,$3,$4,$5,-9}' $GENO_FILE$CHR_MAX".fam" > $FAM_FILE
echo "Correctly formatted FAM file written on: "
echo $FAM_FILE

