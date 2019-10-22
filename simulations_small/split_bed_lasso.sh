#!/bin/bash

################
# Parse input  #
################

RESOLUTION=$1
FOLD=$2

#########
# Setup #
#########

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

echo "----------------------------------------------------------------------------------------------------"
echo "Merging knockoff-augmented BED files the following parameters:"
echo "  - Chromosome range  :" $CHR_MIN".."$CHR_MAX
echo "  - Resolution        :" $RESOLUTION
echo "  - FOLD              :" $FOLD
echo "----------------------------------------------------------------------------------------------------"
echo ""

#################################
# Initialize input/output paths #
#################################

# Input
PHENO_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small/phenotypes"
PHENO_FILE=$PHENO_DIR"/"$EXPERIMENT"_phenotypes_s"$N_SIGNALS"_"$FOLD".tab"
GENO_FILE="/scratch/PI/candes/ukbiobank_tmp/knockoffs/"$RESOLUTION"_K50/ukb_gen_chr"

# List of individuals
FAM_FILE="/scratch/PI/candes/ukbiobank_tmp/simulations_small/subjects/ukb_gen.fam"

###################################
# Split large BED file into folds #
###################################

# Location of BED file for all chromosomes
DATA_BIG_DIR="/scratch/PI/candes/ukbiobank_tmp/augmented_data_big"
DATA_BIG=$DATA_BIG_DIR"/ukb_gen_"$RESOLUTION

# Location of output
LASSO_DATA_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small/lasso_data/"
LASSO_DATA=$LASSO_DATA_DIR"/ukb_gen_"$RESOLUTION

FOLD_INCLUDE_FILE="/scratch/PI/candes/ukbiobank_tmp/simulations_small/subjects/ukb_gen_"$FOLD".fam"

plink \
  --bfile $DATA_BIG \
  --keep $FOLD_INCLUDE_FILE \
  --make-bed \
  --memory 39000 \
  --out $LASSO_DATA"_fold_"$FOLD

rm $LASSO_DATA"_fold_"$FOLD".log"
