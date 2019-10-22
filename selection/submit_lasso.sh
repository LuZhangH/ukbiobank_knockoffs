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

BAD_CPUS="sh-105-22"
RESET=1

# Slurm parameters
PART=candes       # Partition names
MEMO=100G          # Memory required
TIME=00-08:00:00  # Time required
CORE=10           # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_CPUS"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Parameters
CLUMPING_LIST=("Radj2" "Radj5") # "Radj2" "Radj5" "Radj10" "Radj20" "Radj50"
PHENOTYPE_LIST=("rd") # "hypothyroidism" "cvd" "rd"

# Convert chromosomes, one-by-one
for CLUMPING in "${CLUMPING_LIST[@]}"; do
  for PHENOTYPE in "${PHENOTYPE_LIST[@]}"; do
    OUTFILE="/scratch/PI/candes/ukbiobank_tmp/association/"$PHENOTYPE"/"$CLUMPING"_K50/lasso.txt"
    if [[ ! -f $OUTFILE ]] || [[ $RESET -eq 1 ]]; then
      # Script to be run
      SCRIPT="select_lasso.sh $PHENOTYPE $CLUMPING"
      # Define job name for this chromosome
      JOBN="lasso_"$PHENOTYPE"_"$CLUMPING
      OUTF=$LOGS"/"$JOBN".out"
      ERRF=$LOGS"/"$JOBN".err"
      # Assemble slurm order for this job
      ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
      # Print order
      echo $ORD
      # Submit order
      $ORD
    else
      echo "File "$OUTFILE" exists"
    fi
  done
done
