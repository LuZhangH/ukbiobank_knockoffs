#!/usr/bin/env bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Prepare phased haplotype for use with fastPHASE
#
# wrapper script for job submission on Broad UGER cluster
#
# adopted from Nik and Matteo's scripts
#
# Author: Lu Zhang
# Date:   12/30/2019

# Slurm parameters
TASK=conv          # Task description (short string)
#PART=candes        # Partition names
MEMO=200G           # Memory required
TIME=01-00:00:00   # Time required (4d)
CORE=1             # Cores required

# Assemble order prefix
ORDP="qsub -l h_rt="$TIME" -l h_vmem="$MEMO" -pe smp "$CORE" -R y -binding linear:"$CORE"" 
#ORDP="sbatch --mem="$MEMO" -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=~/ukbiobank_knockoffs/fastphase/logs
#LOGS=logs
mkdir -p $LOGS

#CHR_LIST=(1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
#CHR_LIST=(1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
CHR_LIST=(21)
# Convert chromosomes, one-by-one
for CHR in "${CHR_LIST[@]}"; do
  # Run this command
  SCRIPT="convert_v2.sh -c $CHR -n 1000 -b 10000"

  # Define job name for this chromosome
  JOBN=$TASK"_"$CHR
  OUTF=$LOGS"/"$JOBN".out"
  ERRF=$LOGS"/"$JOBN".err"
  # Assemble slurm order for this job
  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
  # Print order
  echo $ORD
  # Submit order
  $ORD
done  
