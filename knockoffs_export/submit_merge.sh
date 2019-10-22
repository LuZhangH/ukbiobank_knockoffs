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
RESOLUTION_LIST=("Radj100" "Radj75" "Radj50" "Radj20" "Radj10" "Radj5" "Radj2")
SEED_LIST=(1)

# Slurm parameters
PART=candes              # Partition names
MEMO=20G                 # Memory required (60GB for merging BED)
TIME=02-00:00:00         # Time required (6h for for merging BED)
CORE=1                   # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
  for SEED in "${SEED_LIST[@]}"; do
    OUT_DIR="/scratch/PI/candes/ukbiobank_tmp/augmented_data_big"
    OUT_FILE=$OUT_DIR"/ukb_gen_"$RESOLUTION"_s"$SEED".bed"
    #if [ ! -s $OUT_FILE ]; then
    # Script to be run
    SCRIPT="merge_bed.sh $RESOLUTION $SEED"
    # Define job name for this chromosome
    JOBN="merge_"$RESOLUTION"_s"$SEED
    OUTF=$LOGS"/"$JOBN".out"
    ERRF=$LOGS"/"$JOBN".err"
    # Assemble slurm order for this job
    ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
    # Print order
    echo $ORD
    # Submit order
    $ORD
    #fi
  done
done
