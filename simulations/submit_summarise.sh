#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Export augmented genetic data to shared folder.
# 
# Author: Matteo Sesia
# Date:   01/15/2019

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Slurm parameters
PART=candes,hns   # Partition names
#MEMO=40G             # Memory required (40 GB for knockoffs)
#TIME=00-1:00:00      # Time required (1h)
#CORE=8              # Cores required (8 cores for knockoffs)
MEMO=5G
TIME=00-00:20:00
CORE=1

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

EXPERIMENT_LIST=("one" "two" "three" "four" "five" "six" "seven" "eight" "nine" "ten")
#EXPERIMENT_LIST=("six" "seven")
SUMMARY="multires"
#SUMMARY="knockoffs"
#SUMMARY="lmm_fine"
#SUMMARY="lmm_locus"

for EXPERIMENT in "${EXPERIMENT_LIST[@]}"; do    

  # Script to be run
  SCRIPT="summarise.sh $SUMMARY $EXPERIMENT 1"
  # Define job name for this chromosome
  JOBN="summary_"$SUMMARY"_"$EXPERIMENT
  OUTF=$LOGS"/"$JOBN".out"
  ERRF=$LOGS"/"$JOBN".err"
  # Assemble slurm order for this job
  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
  # Print order
  echo $ORD
  # Submit order
  $ORD

done
