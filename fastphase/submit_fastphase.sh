#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Run fastPHASE on each chromosome separately.
#
# Author: Matteo Sesia
# Date:   06/20/2018

# Slurm parameters
TASK=FP                 # Task description (short string)
PART=candes             # Partition names
MEMO=50G                # Memory required (50G for chr 1)
TIME=6-00:00:00        # Time required (6d for K=50)
CORE=1                  # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# List of parameters
CHR_LIST=$(seq 1 22)
K_LIST=(100) # 1 2 5 10 15 50 75 100 125 150 175 200

# Convert chromosomes, one-by-one
for K in ${K_LIST[@]}; do
for CHR in $CHR_LIST; do
  # Run this command
  SCRIPT="fastphase.sh -c $CHR -k $K"

  # Define job name for this chromosome
  JOBN=$TASK$CHR"K"$K
  OUTF=$LOGS"/"$JOBN".out"
  ERRF=$LOGS"/"$JOBN".err"
  # Assemble slurm order for this job
  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
  # Print order
  echo $ORD
  # Submit order
  $ORD
done
done
