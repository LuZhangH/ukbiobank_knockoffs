#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Export augmented genetic data to shared folder.
# 
# Author: Matteo Sesia
# Date:   09/09/2018

# Parameters
#CHR_LIST=$(seq 1 22)
#K_LIST=(50)
#CLUMPING_LIST=(100 75 50 20 10 5 2 1 0.5)
CHR_LIST=$(seq 1 22)
K_LIST=(50)
CLUMPING_LIST=(100 75 50 20 10 5 2) #100 
SEED_LIST=(1)

# Convert chromosomes, one-by-one
for SEED in ${SEED_LIST[@]}; do
  for CLUMPING in ${CLUMPING_LIST[@]}; do
    for CHR in $CHR_LIST; do
      for K in ${K_LIST[@]}; do
        # Script to be run
        COMMAND="./export.sh -c $CHR -m Radj -r $CLUMPING -k $K -s $SEED --reset"
        # Print order
        echo $COMMAND
        # Submit order
        $COMMAND
      done
    done
  done
done
