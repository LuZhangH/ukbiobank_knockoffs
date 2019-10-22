#!/bin/bash

################
# Parse input  #
################

RESOLUTION=$1
SEED=$2

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
echo "----------------------------------------------------------------------------------------------------"
echo ""

#############################
# Merge BED files for lasso #
#############################

# Input
GENO_FILE="/scratch/PI/candes/ukbiobank_tmp/knockoffs/"$RESOLUTION"_K50/ukb_gen_chr"

# Location for temporary data storage
LASSO_DATA_DIR="/scratch/PI/candes/ukbiobank_tmp/augmented_data_big"
mkdir -p $LASSO_DATA_DIR

# File name for temporary data storage
LASSO_DATA=$LASSO_DATA_DIR"/ukb_gen_"$RESOLUTION"_s"$SEED

if [ ! -s $LASSO_DATA".bed" ]; then
  # Merge multiple BED files
  GENO_FILE="/scratch/PI/candes/ukbiobank_tmp/knockoffs/"$RESOLUTION"_K50_s"$SEED"/ukb_gen_chr"
  MERGE_LIST=$LASSO_DATA_DIR"/"$RESOLUTION"_mergelist.txt"
  for CHR in $CHR_LIST; do
    if [ $CHR == 1 ]; then
      rm -rf $MERGE_LIST
      touch $MERGE_LIST
    else
      echo $GENO_FILE$CHR >> $MERGE_LIST
    fi
  done

  plink \
    --bfile $GENO_FILE"1" \
    --merge-list $MERGE_LIST \
    --make-bed \
    --memory 39000 \
    --out $LASSO_DATA

  rm $LASSO_DATA".log"
else
  echo "--------------------------------------------------"
  echo "Skipping merge because:"
  echo " "$LASSO_DATA".bed exists."
  echo "--------------------------------------------------"
fi

##################################
# Make list of original variants #
##################################

Rscript --vanilla list_original.R $RESOLUTION

######################
# Convert BED to FBM #
######################

# Call R script to convert BED to FBM
if [ ! -s $LASSO_DATA".bk" ] | [ ! -s $LASSO_DATA".rds" ]; then
  echo "--------------------------------------------------"
  echo "Converting BED to FBM"
  echo "--------------------------------------------------"
  source activate ukb
  Rscript --vanilla make_FBM.R $RESOLUTION $SEED
else
  echo "--------------------------------------------------"
  echo "Skipping conversion of BED to FBM because:"
  echo " "$LASSO_DATA".bk exists."
  echo " "$LASSO_DATA".rds exists."
  echo "--------------------------------------------------"
fi
