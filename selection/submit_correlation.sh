#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Compute correlation between variables and knockoffs.
#
# Author: Matteo Sesia
# Date:   10/01/2018

# Parameters
CHR_LIST=$(seq 1 22)
CLUMPING_LIST=("Radj50" "Radj20" "Radj10" "Radj5" "Radj2")
#CHR_LIST=$(seq 22 22)
#CLUMPING_LIST=("Radj20")

# Slurm parameters
PART=candes       # Partition names
MEMO=50G           # Memory required
TIME=00-00:20:00  # Time required
CORE=1            # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Convert chromosomes, one-by-one
# Convert chromosomes, one-by-one
for CLUMPING in "${CLUMPING_LIST[@]}"; do
    for CHR in $CHR_LIST; do
      # Script to be run
      SCRIPT="compute_correlation.sh $CHR $CLUMPING"
      # Define job name for this chromosome
      JOBN="chr"$CHR"_"$CLUMPING
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
