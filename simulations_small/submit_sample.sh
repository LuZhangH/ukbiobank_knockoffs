#!/bin/bash

BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
BAD_NODES=$(cat $BAD_NODES_FILE)
echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Slurm parameters
PART=candes,hns,owners        # Partition names
MEMO=100G           # Memory required
TIME=00-02:00:00   # Time required
CORE=1            # Cores required

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Script to be run
SCRIPT="sample.sh"
# Define job name for this chromosome
JOBN="sample"
OUTF=$LOGS"/"$JOBN".out"
ERRF=$LOGS"/"$JOBN".err"
# Assemble slurm order for this job
ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
# Print order
echo $ORD
# Submit order
$ORD
