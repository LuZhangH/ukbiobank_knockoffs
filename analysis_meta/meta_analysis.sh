#!/bin/bash
#
# Cross-reference with knockoffs
#

# Choose phenotype
PHENOTYPE=$1 # body_HEIGHTz
CLUMP_THRESHOLD=$2 # 0.00000005
USE_SMALL=$3 # 0

mkdir -p "/scratch/PI/candes/ukbiobank_tmp/meta/summary"

source activate ukb
#Rscript --vanilla meta_analysis.R $PHENOTYPE $CLUMP_THRESHOLD $USE_SMALL
Rscript --vanilla meta_analysis_stable.R $PHENOTYPE $CLUMP_THRESHOLD $USE_SMALL
