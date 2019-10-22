#!/bin/bash

################
# Parse input  #
################

#########
# Setup #
#########

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

echo "----------------------------------------------------------------------------------------------------"
echo "Computing PCA:"
echo "  - Chromosome range  :" $CHR_MIN".."$CHR_MAX
echo "----------------------------------------------------------------------------------------------------"
echo ""

#################################
# Initialize input/output paths #
#################################

GEN_DIR="/scratch/groups/candes/ukbiobank_tmp/augmented_data_big"
GENO_FILE=$GEN_DIR"/ukb_gen_Radj100"
SNP_EXCLUDE=$GEN_DIR"/ukb_gen_Radj100_knockoff.txt"

# Output dir
PCA_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/pca"
mkdir -p $PCA_DIR
PCA_FILE=$PCA_DIR"/ukb_gen"

plink --bfile $GENO_FILE \
      --exclude $SNP_EXCLUDE \
      --indep-pairwise 50 100 0.5 \
      --pca \
      --out $PCA_FILE

