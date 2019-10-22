#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Merge augmented genetic data
# 
# Author: Matteo Sesia
# Date:   02/08/2018

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Parameters
RESOLUTION_LIST=("Radj100" "Radj75" "Radj50" "Radj20" "Radj10" "Radj5" "Radj2" "Radj1")
#RESOLUTION_LIST=("Radj1")
FOLD_LIST=$(seq -f "%02g" 1 10)

# Slurm parameters
PART=candes,owners       # Partition names
MEMO=40G                 # Memory required
TIME=00-6:00:00          # Time required
CORE=1                   # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
  for FOLD in $FOLD_LIST; do
    OUT_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small/lasso_data"
    OUT_FILE=$OUT_DIR"/ukb_gen_"$RESOLUTION"_fold_"$FOLD".bed"
    if [ ! -s $OUT_FILE ]; then
      # Script to be run
      SCRIPT="split_bed_lasso.sh $RESOLUTION $FOLD"
      # Define job name for this chromosome
      JOBN="split_lasso_"$RESOLUTION"_"$FOLD
      OUTF=$LOGS"/"$JOBN".out"
      ERRF=$LOGS"/"$JOBN".err"
      # Assemble slurm order for this job
      ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
      # Print order
      echo $ORD
      # Submit order
      #$ORD
    fi
  done
done
