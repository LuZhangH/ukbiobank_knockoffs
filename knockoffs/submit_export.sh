#!/bin/bash
# UK Biobank GWAS
#
# Class: SLURM dispatcher
#
# Export knockoff-augmented genotype data.
#
# Author: Matteo Sesia
# Date:   06/20/2018

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Slurm parameters
TASK=export           # Task description
PART=candes      # Partition names
MEMO=5G          # Memory required
TIME=00-12:00:00 # Time required
CORE=1           # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Script to be run
SCRIPT="export_all.sh"

# Define job name for this chromosome
JOBN=$TASK
OUTF=$LOGS"/"$JOBN".out"
ERRF=$LOGS"/"$JOBN".err"
# Assemble slurm order for this job
ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
# Print order
echo $ORD
# Submit order
$ORD
