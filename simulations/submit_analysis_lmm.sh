#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Export augmented genetic data to shared folder.
# 
# Author: Matteo Sesia
# Date:   09/09/2018
#
# Architecture:
#  /scratch/groups/candes/ukbiobank_tmp/simulations/phenotypes/architecture.txt
#

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Parameters
FOLD_LIST=$(seq -f "%g" 1 1)

# Slurm parameters
PART=candes,hns,owners,normal          # Partition names
MEMO=150G                 # Memory required (50G for BOLT, 100GB for susie)
TIME=01-00:00:00         # Time required (2d for BOLT)
CORE=1                   # Cores required (8 for BOLT)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"
ORDP=$ORDP" --mail-user=matteo1992@gmail.com --mail-type=END"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Loop over configurations and chromosomes
cat parameters_h2.txt | while read PARAM_LINE
do 
  PARAM_LINE=($PARAM_LINE)
  EXPERIMENT_NAME=${PARAM_LINE[0]}
  EXPERIMENT_NUMB=${PARAM_LINE[1]}
  for FOLD in $FOLD_LIST; do
    echo "Preparing experiment: ("$EXPERIMENT_NAME","$EXPERIMENT_NUMB"), fold: "$FOLD

    #OUT_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/bolt/fold_"$FOLD
    #OUT_FILE=$OUT_DIR"/stats_"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB".txt"
    OUT_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/susie/fold_"$FOLD
    OUT_FILE=$OUT_DIR"/"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB"_0.00000005_0.9_nobias.txt"
    #OUT_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/caviar/fold_"$FOLD
    #OUT_FILE=$OUT_DIR"/"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB"_0.00000005.txt"
    #OUT_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/BH/fold_"$FOLD
    #OUT_FILE=$OUT_DIR"/selections_"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB".txt1"
    if [ ! -s $OUT_FILE ]; then
      # Script to be run
      SCRIPT="analysis_lmm.sh $EXPERIMENT_NAME $EXPERIMENT_NUMB $FOLD"
      # Define job name for this chromosome
      JOBN="lmm_"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB"_"$FOLD
      OUTF=$LOGS"/"$JOBN".out"
      ERRF=$LOGS"/"$JOBN".err"
      # Assemble slurm order for this job
      ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
      # Print order
      echo $ORD
      # Submit order
      #$ORD
      # Wait between submissions because BOLT-LMM messes up if starting at the same time
      #sleep 1      
    fi
  done
done
