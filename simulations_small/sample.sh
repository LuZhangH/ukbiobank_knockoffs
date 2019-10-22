#!/bin/bash

# Parameters
FOLD_SIZE=30000

# Create output directories
mkdir -p "/scratch/groups/candes/ukbiobank_tmp/simulations_small/phenotypes/"

# Sample Y|X
source activate ukb
#Rscript --vanilla sample.R

# Split phenotype files into folds
#EXPERIMENT_LIST=("one" "two" "three" "four" "five" "six" "seven" "eight" "nine" "ten")
EXPERIMENT_LIST=("sevenps")
for EXPERIMENT in "${EXPERIMENT_LIST[@]}"; do    

  TMP_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small"
  OUT_DIR=$TMP_DIR"/phenotypes"
  PHENO_FILE=$TMP_DIR"/phenotypes/"$EXPERIMENT"_phenotypes.tab"

  echo "Splitting phenotype file $PHENO_FILE into folds of size $FOLD_SIZE ..."
  split <(tail -n +2 $PHENO_FILE) -l $FOLD_SIZE --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".tab" $TMP_DIR"/phenotypes/"$EXPERIMENT"_phenotypes_"

  HEADER=$(head -n 1 $PHENO_FILE)
  for FILE in $TMP_DIR"/phenotypes/"$EXPERIMENT"_phenotypes_"*.tab; do
    echo "Written file $FILE"
    echo -e "$HEADER\n$(cat $FILE)" > $FILE
  done

done

# Split individuals into folds
#./make_folds.sh $FOLD_SIZE
#./bolt_init.sh
