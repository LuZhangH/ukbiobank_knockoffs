#!/bin/bash
#
# Clump Po-Ru's p-values
#

# Choose phenotype
PHENOTYPE=$1 # body_HEIGHTz
CLUMP_THRESHOLD=$2 # 0.00000005

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

# Main storage location
META_DIR="/scratch/groups/candes/ukbiobank_tmp/meta"
mkdir -p $META_DIR
mkdir -p $META_DIR"/stats/"

# Storage location
CLUMP_DIR=$META_DIR"/clumped"
mkdir -p $CLUMP_DIR

# Define genetic files
GENO_FILE="/scratch/groups/candes/ukbiobank/genotypes/ukb_gen_chr"

# Define stats file
STATS_FILE=$META_DIR"/stats/"$PHENOTYPE".sumstats"

for CHR in $CHR_LIST; do
  # Clumping of LMM p-values
  CLUMP_FILE=$CLUMP_DIR"/tmp/"$PHENOTYPE"_"$CLUMP_THRESHOLD
  if [ ! -f $CLUMP_FILE"_chr"$CHR".clumped" ]; then
    echo "Clumping "$PHENOTYPE" (chromosome $CHR at threshold "$CLUMP_THRESHOLD")"
    echo $CLUMP_FILE"_chr"$CHR".clumped"
    touch $CLUMP_FILE"_chr"$CHR".clumped"
    plink --bfile $GENO_FILE$CHR \
          --clump $STATS_FILE \
          --clump-p1 $CLUMP_THRESHOLD \
          --clump-r2 0.01 \
          --clump-kb 5000 \
          --out $CLUMP_FILE"_chr"$CHR
    rm $CLUMP_FILE"_chr"$CHR".log"
    rm $CLUMP_FILE"_chr"$CHR".nosex"
  else
    echo "Skipping clumping of "$PHENOTYPE" (chromosome $CHR at threshold "$CLUMP_THRESHOLD")"
  fi
done

# Combine results of clumping into a single file
echo "Combining clumping files for $PHENOTYPE at threshold $CLUMP_THRESHOLD"
CLUMP_FILE=$CLUMP_DIR"/"$PHENOTYPE"_"$CLUMP_THRESHOLD".clumped"
for CHR in $CHR_LIST; do
  CLUMP_CHR_FILE=$CLUMP_DIR"/tmp/"$PHENOTYPE"_"$CLUMP_THRESHOLD"_chr"$CHR".clumped"
  if [[ $CHR == 1 ]]; then
    cat $CLUMP_CHR_FILE > $CLUMP_FILE
  else
    tail -n +2 $CLUMP_CHR_FILE >> $CLUMP_FILE  
  fi
  echo "Appended clumps from chromosome "$CHR
done
# Remove empty lines
sed -i '/^$/d' $CLUMP_FILE
echo "Results written on "$CLUMP_FILE
