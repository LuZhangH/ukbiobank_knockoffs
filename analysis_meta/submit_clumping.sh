#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Clump the p-values computed by Po-Ru
#
# Author: Matteo Sesia
# Date:   29/03/2019

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Parameters
PHENOTYPE_LIST=("body_HEIGHTz" "body_BMIz" "blood_PLATELET_COUNT" "bp_SYSTOLICadjMEDz" "disease_CARDIOVASCULAR" "disease_HYPOTHYROIDISM_SELF_REP")
#PHENOTYPE_LIST=("body_HEIGHTz")

CLUMP_THRESHOLD_LIST=("0.00000005") # "0.000000005"

# Slurm parameters
PART=candes,hns,normal,owners       # Partition names
MEMO=5G              # Memory required
TIME=00-01:00:00      # Time required
CORE=1                # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Loop over configurations and chromosomes
for CLUMP_THRESHOLD in "${CLUMP_THRESHOLD_LIST[@]}"; do
  for PHENOTYPE in "${PHENOTYPE_LIST[@]}"; do
    for USE_SMALL in $(seq 0 0); do # This is not for clumping
      for USE_BH in $(seq 0 0); do # This is not for clumping
        #OUT_DIR="/scratch/groups/candes/ukbiobank_tmp/meta/summary"
        #OUT_FILE=$OUT_DIR"/"$PHENOTYPE"_stats.txt"
        # if [ ! -s $OUT_FILE ]; then
        # Script to be run
        #SCRIPT="clump.sh $PHENOTYPE $CLUMP_THRESHOLD"
        SCRIPT="meta_analysis.sh $PHENOTYPE $USE_SMALL $USE_BH"
        # Define job name for this chromosome
        JOBN="clump_"$PHENOTYPE"_"$CLUMP_THRESHOLD"_"$USE_SMALL
        OUTF=$LOGS"/"$JOBN".out"
        ERRF=$LOGS"/"$JOBN".err"
        # Assemble slurm order for this job
        ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
        # Print order
        echo $ORD
        # Submit order
        $ORD
        # else
        #   echo "Skipping analysis of "$PHENOTYPE" because"
        #   echo $OUT_FILE" exists and is not empty"
        # fi
      done
    done
  done
done
