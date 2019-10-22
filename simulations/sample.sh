#!/bin/bash

# Input parameters
EXPERIMENT=$1
N_SIGNALS=$2

# Create output directories
mkdir -p "/scratch/groups/candes/ukbiobank_tmp/simulations/phenotypes/"

# Sample Y|X
source activate ukb
Rscript --vanilla sample.R

# # Copy family file with individiduals that passed QC
# FAM_ORIGINAL="/scratch/groups/candes/ukbiobank_tmp/knockoffs/Radj100_K50/ukb_gen_chr1.fam"
# FAM_SIMULATIONS="/scratch/groups/candes/ukbiobank_tmp/simulations/subjects/ukb_gen_00.fam"
# cp $FAM_ORIGINAL $FAM_SIMULATIONS

# # Initialize auxiliary files for BOLT-LMM
# ./bolt_init.sh
