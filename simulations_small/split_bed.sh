#!/bin/bash

################
# Parse input  #
################

CHR=$1

#########
# Setup #
#########

# Range of chromosomes to include in the analysis
FOLD_LIST=$(seq -f "%02g" 1 10)

echo "----------------------------------------------------------------------------------------------------"
echo "Splitting genotype BED files the following parameters:"
echo "  - Chromosome    :" $CHR
echo "----------------------------------------------------------------------------------------------------"
echo ""

#################################
# Initialize input/output paths #
#################################

# Input
GENO_FILE_IN="/scratch/groups/candes/ukbiobank/genotypes/ukb_gen_chr"

# Output
GENO_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations_small/bolt_data"
mkdir -p $GENO_DIR
GENO_FILE_OUT=$GENO_DIR"/ukb_gen"

# List of subjects
FAM_FILE="/scratch/groups/candes/ukbiobank_tmp/simulations_small/subjects/ukb_gen.fam"

# List of variants
SNP_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations_small/variants"
SNP_EXCLUDE_FILE=$SNP_DIR"/ukb_gen_exclude.txt"

###################################
# Split large BED file into folds #
###################################

for FOLD in $FOLD_LIST; do

  FOLD_FILE="/scratch/PI/candes/ukbiobank_tmp/simulations_small/subjects/ukb_gen_"$FOLD".fam"
  REMOVE_FILE="/scratch/PI/candes/ukbiobank_tmp/simulations_small/subjects/ukb_gen_"$FOLD"_exclude.fam"

  plink \
      --bfile $GENO_FILE_IN$CHR \
      --keep $FOLD_FILE \
      --exclude $SNP_EXCLUDE_FILE \
      --make-bed \
      --out $GENO_FILE_OUT"_fold"$FOLD"_chr"$CHR

  rm $GENO_FILE_OUT"_fold"$FOLD"_chr"$CHR".log" $GENO_FILE_OUT"_fold"$FOLD"_chr"$CHR".nosex"

done
