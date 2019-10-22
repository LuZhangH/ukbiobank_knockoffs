#!/bin/bash

FOLD_SIZE=$1

# Output directory
DAT_DIR="/scratch/PI/candes/ukbiobank_tmp"
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small"
mkdir -p $TMP_DIR
mkdir -p $TMP_DIR"/subjects"

# Split list of individuals into folds
FAM_FILE=$DAT_DIR"/augmented_data/Radj10_K50/ukb_gen_chr1.fam"
NUM_SAMPLES=$(wc -l < "$FAM_FILE")
echo "Total number of individuals: "$NUM_SAMPLES
FOLD_NUM=$(($NUM_SAMPLES / $FOLD_SIZE))
echo "Splitting into "$FOLD_NUM" folds"
split $FAM_FILE -l $FOLD_SIZE --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".fam" $TMP_DIR"/subjects/ukb_gen_"
echo "Saved lists of individuals in: "$TMP_DIR"/subjects/"
