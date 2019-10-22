#!/bin/bash

# Input parameters
CHR=$1
CLUMPING=$2

# Fixed parameters
K=50

# Input and output files
BASENAME="/scratch/PI/candes/ukbiobank_tmp/augmented_data/"$CLUMPING"_K"$K"/ukb_gen_chr"$CHR
OUT_BASENAME="/scratch/PI/candes/ukbiobank_tmp/augmented_data/"$CLUMPING"_K"$K"/ukb_gen_chr"$CHR

# Load plink and datamash
export PATH=$PI_HOME/bin/:$PATH

plink --bfile $BASENAME  \
      --r --ld-window 6 \
      --out $OUT_BASENAME

# Remove log file
rm $OUT_BASENAME".log"
