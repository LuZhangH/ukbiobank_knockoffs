#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Run fastPHASE on masked data to tune K.
#
# Author: Matteo Sesia
# Date:   07/17/2018

# Slurm parameters
TASK=CV            # Task description (short string)
PART=candes,hns,owners        # Partition names
MEMO=10G           # Memory required
TIME=00-06:00:00   # Time required
CORE=1             # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# List of parameters
CHR_LIST=$(seq 22 22)
K_LIST=(1 2 5 10 20 30 50 75 100 150 200)
MASK=01
N=100000

# Convert chromosomes, one-by-one
for K in ${K_LIST[@]}; do
for CHR in $CHR_LIST; do
  # Run this command
  SCRIPT="CV_fastphase.sh -c $CHR -k $K -n $N -r"

  # Define job name for this chromosome
  JOBN=$TASK$CHR"K"$K"N"$N
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
