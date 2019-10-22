#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Download the p-values computed by Po-Ru
# 
# Author: Matteo Sesia
# Date:   29/03/2019

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Slurm parameters
PART=candes,hns       # Partition names
MEMO=10G              # Memory required
TIME=00-04:00:00      # Time required
CORE=1                # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Script to be run
SCRIPT="download.sh"
# Define job name for this chromosome
JOBN="download"
OUTF=$LOGS"/"$JOBN".out"
ERRF=$LOGS"/"$JOBN".err"
# Assemble slurm order for this job
ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
# Print order
echo $ORD
# Submit order
$ORD
