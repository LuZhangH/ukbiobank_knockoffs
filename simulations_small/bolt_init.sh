#!/bin/bash

# List of folds
FOLD_LIST=$(seq -f "%02g" 1 10)

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22

echo "Initializing list of variants and individuals for BOLT-LMM on $FOLD"
echo "  First CHR         :" $CHR_MIN
echo "  Last CHR          :" $CHR_MAX

###############################
# Make list of SNPs to exlude #
###############################
# Input
GENO_FILE="/scratch/PI/candes/ukbiobank/genotypes/ukb_gen_chr"

# Make directory for output
SNP_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small/variants"
mkdir -p $SNP_DIR

# Work on one chromosome at a time
for CHR in $(seq $CHR_MIN $CHR_MAX); do
  echo "Making list of excluded SNPs from chromosome $CHR ..."
  # Remove knockoffs from the list of original variables used in this experiment
  BIM_FILE="/scratch/PI/candes/ukbiobank_tmp/augmented_data/Radj50_K50/ukb_gen_chr"$CHR".bim"
  SNP_FILE=$SNP_DIR"/ukb_gen_chr"$CHR".bim"
  awk -F"[\t.]" '{print $2}' $BIM_FILE | uniq > $SNP_FILE
  # List the names of the SNPs in the original data
  SNP_TMP_FILE=$SNP_DIR"/ukb_gen_chr"$CHR"_tmp.txt"
  awk -F"[\t.]" '{print $2}' $GENO_FILE$CHR".bim" > $SNP_TMP_FILE
  # Find the names of the SNPs to exclude
  SNP_EXCLUDE_FILE=$SNP_DIR"/ukb_gen_chr"$CHR"_exclude.txt"
  comm -23 <(sort -i $SNP_TMP_FILE) <(sort -i $SNP_FILE) > $SNP_EXCLUDE_FILE
  # Remove temporary files
  rm $SNP_FILE $SNP_TMP_FILE
  echo "List written on: "$SNP_EXCLUDE_FILE
done

# Combine list of files for different chromosomes
SNP_EXCLUDE_FILE=$SNP_DIR"/ukb_gen_exclude.txt"
cat $SNP_DIR"/ukb_gen_chr"*"_exclude.txt" > $SNP_EXCLUDE_FILE
echo "Final list of exluded variants written on: "
echo $SNP_EXCLUDE_FILE

######################################
# Make list of individuals to exlude #
######################################
for FOLD in $FOLD_LIST; do    
  # List of individuals
  FOLD_FILE="/scratch/PI/candes/ukbiobank_tmp/simulations_small/subjects/ukb_gen_"$FOLD".fam"
  FAM_FILE="/scratch/PI/candes/ukbiobank_tmp/simulations_small/subjects/ukb_gen.fam"
  awk '{print $1,$2,$3,$4,$5,-9}' $GENO_FILE$CHR_MAX".fam" > $FAM_FILE

  # List of excluded individuals
  FOLD_EXCLUDE_FILE="/scratch/PI/candes/ukbiobank_tmp/simulations_small/subjects/ukb_gen_"$FOLD"_exclude.fam"
  comm -23 <(sort -i $FAM_FILE) <(sort -i $FOLD_FILE) > $FOLD_EXCLUDE_FILE
  echo "List of individuals exluded from fold "$FOLD" written on:"
  echo $FOLD_EXCLUDE_FILE
done
