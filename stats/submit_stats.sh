#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Compute variants stats for each chromosome separately.
#
# Author: Matteo Sesia
# Date:   08/13/2018

# Slurm parameters
TASK=stats              # Task description (short string)
PART=candes             # Partition names
MEMO=20G                # Memory required (50G for chr 1)
TIME=00-00:30:00        # Time required
CORE=1                  # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# List of parameters
CHR_LIST=$(seq 1 22)

# Convert chromosomes, one-by-one
for CHR in $CHR_LIST; do
  # Run this command
  SCRIPT="variant_stats.sh -c $CHR"

  # Define job name for this chromosome
  JOBN=$TASK$CHR
  OUTF=$LOGS"/"$JOBN".out"
  ERRF=$LOGS"/"$JOBN".err"
  # Assemble slurm order for this job
  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
  # Print order
  echo $ORD
  # Submit order
  $ORD
done
