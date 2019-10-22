#!/bin/bash

FOLD_SIZE=30000

# Output directory
DAT_DIR="/scratch/PI/candes/ukbiobank_tmp"
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp/analysis"
mkdir -p $TMP_DIR
mkdir -p $TMP_DIR"/individuals_split"

# Split list of individuals into folds
FAM_FILE=$DAT_DIR"/augmented_data/Radj10_K50/ukb_gen_chr1.fam"
NUM_SAMPLES=$(wc -l < "$FAM_FILE")
echo "Total number of individuals: "$NUM_SAMPLES
FOLD_NUM=$(($NUM_SAMPLES / $FOLD_SIZE))
echo "Splitting into "$FOLD_NUM" folds"
split $FAM_FILE -l $FOLD_SIZE --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".fam" $TMP_DIR"/individuals_split/ukb_gen_"
echo "Saved lists of individuals in: "$TMP_DIR"/individuals_split/"

######################################
# Make list of individuals to exlude #
######################################
FOLD_LIST=$(seq -f "%02g" 1 1)

for FOLD in $FOLD_LIST; do    
  # List of individuals
  FOLD_FILE="/scratch/PI/candes/ukbiobank_tmp/analysis/individuals_split/ukb_gen_"$FOLD".fam"  
  FAM_FILE="/scratch/groups/candes/ukbiobank_tmp/analysis/individuals/ukb_individuals.fam"
  FAM_NEW_FILE="/scratch/PI/candes/ukbiobank_tmp/analysis/individuals_split/ukb_gen.fam"
  awk '{print $1,$2,$3,$4,$5,-9}' $FAM_FILE > $FAM_NEW_FILE

  # List of excluded individuals
  FOLD_EXCLUDE_FILE="/scratch/PI/candes/ukbiobank_tmp/analysis/individuals_split/ukb_gen_"$FOLD"_exclude.fam"
  comm -23 <(sort -i $FAM_NEW_FILE) <(sort -i $FOLD_FILE) > $FOLD_EXCLUDE_FILE
  echo "List of individuals exluded from fold "$FOLD" written on:"
  echo $FOLD_EXCLUDE_FILE
done
