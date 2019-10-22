#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Export augmented genetic data to shared folder.
#
# Author: Matteo Sesia
# Date:   09/09/2018

BAD_CPUS="sh-105-22"

# Parameters
CHR_LIST=$(seq 1 22)
CLUMPING_LIST=("Radj5" "Radj10") # "Radj5" "Radj10" "Radj20" "Radj50") # "Radj2" "Radj5" "Radj10" "Radj20" "Radj50"
PHENOTYPE_LIST=("rd") # "hypothyroidism" "cvd" "rd"

# Slurm parameters
PART=candes       # Partition names
MEMO=50G           # Memory required
TIME=01-00:00:00  # Time required
CORE=1            # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_CPUS"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Convert chromosomes, one-by-one
# Convert chromosomes, one-by-one
for CLUMPING in "${CLUMPING_LIST[@]}"; do
  for PHENOTYPE in "${PHENOTYPE_LIST[@]}"; do
    for CHR in $CHR_LIST; do
      # Script to be run
      SCRIPT="plink_lasso.sh -c $CHR -m $CLUMPING -k 50 -p $PHENOTYPE"
      # Define job name for this chromosome
      JOBN=$PHENOTYPE"_chr"$CHR"_"$CLUMPING
      OUTF=$LOGS"/"$JOBN".out"
      ERRF=$LOGS"/"$JOBN".err"
      # Assemble slurm order for this job
      ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
        # Print order
      echo $ORD
      # Submit order
      #$ORD
    done
  done
done
