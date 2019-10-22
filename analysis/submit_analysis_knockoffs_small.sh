#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Run BOLT-LLM on knockoff-augmented data
#
# Author: Matteo Sesia
# Date:   01/25/2019

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Parameters
PHENOTYPE_LIST=("height" "bmi" "sbp" "platelet" "cvd" "diabetes" "hypothyroidism" "respiratory" "glaucoma")
#PHENOTYPE_LIST=("height")
#RESOLUTION_LIST=("Radj2" "Radj5" "Radj10" "Radj20" "Radj50" "Radj75" "Radj100")
RESOLUTION_LIST=("Radj2")

# Slurm parameters
PART=candes           # Partition names
MEMO=20G                  # Memory required (20GB)
TIME=01-00:00:00          # Time required (8h)
CORE=10                   # Cores required (10)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

# Loop over configurations and chromosomes
for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
  for PHENOTYPE in "${PHENOTYPE_LIST[@]}"; do
    OUT_DIR="/scratch/PI/candes/ukbiobank_tmp/analysis/knockoffs"
    OUT_FILE=$OUT_DIR"/"$PHENOTYPE"_"$RESOLUTION"_fold_01_lasso.txt1"
    
    if [ -f $OUT_FILE ]; then
      echo $OUT_FILE" was found"
    else
      # Script to be run
      SCRIPT="analysis_knockoffs_small.sh $PHENOTYPE $RESOLUTION"
      # Define job name for this chromosome
      JOBN="sknockoffs_"$PHENOTYPE"_"$RESOLUTION
      OUTF=$LOGS"/"$JOBN".out"
      ERRF=$LOGS"/"$JOBN".err"
      # Assemble slurm order for this job
      ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
      # Print order
      echo $ORD
      # Submit order
      $ORD
    fi
  done
done
