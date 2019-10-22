#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Export augmented genetic data to shared folder.
# 
# Author: Matteo Sesia
# Date:   09/09/2018

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Parameters
PHENOTYPE_LIST=("height" "bmi" "sbp" "platelet" "cvd" "diabetes" "hypothyroidism" "respiratory" "glaucoma")
#PHENOTYPE_LIST=("height")

# Slurm parameters
PART=candes,hns       # Partition names
MEMO=80G              # Memory required (80G for BOLT)
TIME=00-12:00:00       # Time required (4d for BOLT)
CORE=10               # Cores required (10 for BOLT)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Loop over configurations and chromosomes
for PHENOTYPE in "${PHENOTYPE_LIST[@]}"; do    
  OUT_DIR="/scratch/groups/candes/ukbiobank_tmp/analysis/bolt"
  OUT_FILE=$OUT_DIR"/"$PHENOTYPE"_stats_fold_01.txt"
  if [ ! -s $OUT_FILE ]; then
    # Script to be run
    SCRIPT="analysis_lmm_small.sh $PHENOTYPE"
    # Define job name for this chromosome
    JOBN="slmm_"$PHENOTYPE
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
  else
    echo "Skipping analysis of "$PHENOTYPE" because"
    echo $OUT_FILE" exists and is not empty"
  fi
done
