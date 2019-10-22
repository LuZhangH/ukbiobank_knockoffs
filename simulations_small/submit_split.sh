#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Split BED files
# 
# Author: Matteo Sesia
# Date:   03/18/2019
#

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Parameters
CHR_LIST=$(seq -f "%g" 1 22)

# Slurm parameters
PART=candes,hns,owners   # Partition names
MEMO=10G                 # Memory required
TIME=00-02:00:00         # Time required
CORE=1                   # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"
ORDP=$ORDP" --mail-user=matteo1992@gmail.com --mail-type=END"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Loop over configurations and chromosomes
for CHR in $CHR_LIST; do
  echo "Preparing split for chromosome "$CHR

  # Script to be run
  SCRIPT="split_bed.sh $CHR"
  # Define job name for this chromosome
  JOBN="split_"$CHR
  OUTF=$LOGS"/"$JOBN".out"
  ERRF=$LOGS"/"$JOBN".err"
  # Assemble slurm order for this job
  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
  # Print order
  echo $ORD
  # Submit order
  $ORD
done
