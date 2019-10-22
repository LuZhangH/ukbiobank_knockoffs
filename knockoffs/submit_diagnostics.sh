#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Generate correlation matrix for one chromosome, using samples
# and variants that passed quality control.
# 
# Generate knockoffs.
#
# Author: Matteo Sesia
# Date:   06/20/2018

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Parameters
CHR_LIST=$(seq 1 21)
K_LIST=(50) #1 2 5 10 15 20 30 50 75 100 125 150 175
CLUMPING_LIST=(100 75 50 20 10 5 2) # 75 50 20 10 5 2 1
SEED_LIST=(0)
#CLUMPING_LIST=(100 75)

# Slurm parameters
TASK=d           # Task description
PART=candes      # Partition names
MEMO=10G         # Memory required (clumping: 200GB, knockoffs: 30GB)
TIME=00-01:00:00 # Time required   (clumping: 4h,  knockoffs: 4d)
CORE=1           # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Convert chromosomes, one-by-one
for SEED in ${SEED_LIST[@]}; do
  for CLUMPING in ${CLUMPING_LIST[@]}; do
    for CHR in $CHR_LIST; do
      for K in ${K_LIST[@]}; do
        OUT_DIR="/scratch/PI/candes/ukbiobank_tmp/knockoff_diagnostics/Radj"$CLUMPING"_K"$K"_s"$SEED
        OUT_FILE=$OUT_DIR"/ukb_gen_chr"$CHR".ld1"
        
        if [ ! -s $OUT_FILE ]; then
          # Script to be run
          SCRIPT="diagnostics.sh -c $CHR -m Radj -r $CLUMPING -k $K -s $SEED"

          # Define job name for this chromosome
          JOBN=$TASK$CHR"K"$K"_"$CLUMPING"_s"$SEED
          OUTF=$LOGS"/"$JOBN".out"
          ERRF=$LOGS"/"$JOBN".err"
          # Assemble slurm order for this job
          ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT

          # Print order
          echo $ORD
          # Submit order
          $ORD
        fi
      done
    done
  done
done
