#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Creates smaller dataset with fewer samples.
# 
# Author: Matteo Sesia
# Date:   02/12/2018

# Input
CHR=$1       # Chromosome number
N_KEEP=$2    # Number of samples to keep

# Load plink
export PATH=$PI_HOME/bin/:$PATH

# Input storage location
DAT_DIR=/scratch/PI/candes/biobank
GEN_DIR=$DAT_DIR/genotypes/EGAD00010001226/001
FAM_DIR=$DAT_DIR/genotypes/sample_app1372

# Ouput storage location
TMP_DIR=/scratch/PI/candes/biobank_tmp
OUT_DIR=$TMP_DIR/data_small
FP_DIR=$TMP_DIR/fastphase_small
mkdir -p $OUT_DIR

# Assemble input file names for this chromosome
BED_FILE=$GEN_DIR"/ukb_cal_chr"$CHR"_v2.bed"
BIM_FILE=$GEN_DIR"/ukb_snp_chr"$CHR"_v2.bim"
FAM_FILE=$FAM_DIR"/ukb1372_cal_chr"$CHR"_v2_s488374.fam"

# Make list of samples to keep
FAM_KEEP=$OUT_DIR/"chr"$CHR"_keep_"$N_KEEP".fam"
awk -v N_KEEP=$N_KEEP 'NR<=N_KEEP {print $1 " " $2}' $FAM_FILE > $FAM_KEEP

# Assemble output file names for this chromosome
BED_NEW=$OUT_DIR"/chr"$CHR"_"$N_KEEP".bed"
BIM_NEW=$OUT_DIR"/chr"$CHR"_"$N_KEEP".bim"
FAM_NEW=$OUT_DIR"/chr"$CHR"_"$N_KEEP".fam"

# Shrink genotypes and save them into the temporary folder
echo "Shrinking chromosome "$CHR" ..."
plink --bed $BED_FILE --bim $BIM_FILE --fam $FAM_FILE \
      --keep $FAM_KEEP \
      --fill-missing-a2 \
      --make-bed --out $OUT_DIR"/chr"$CHR"_"$N_KEEP

# Shrink genotypes and convert them into FASTPHASE format
echo "Shrinking chromosome "$CHR" ..."
plink --bed $BED_FILE --bim $BIM_FILE --fam $FAM_FILE \
      --keep $FAM_KEEP \
      --recode 01 fastphase \
      --out $FP_DIR"/chr"$CHR"_"$N_KEEP".inp"

# Rename FASTPHASE output file
mv $FP_DIR"/chr"$CHR"_"$N_KEEP".inp.chr-"$CHR".recode.phase.inp" $FP_DIR"/chr"$CHR"_"$N_KEEP".inp"

# Remove list of samples to keep
rm $FAM_KEEP
