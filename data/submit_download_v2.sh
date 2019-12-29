#!/usr/bin/env bash

# UK Biobank GWAS
# wrapper script for job submission on Broad UGER cluster
#
# adopted from Nik and Matteo's scripts
#
# Author: Lu Zhang
# Date:   12/28/2019

# BAD_NODES_FILE="/home/users/msesia/sherlock/bad_nodes.txt"
# BAD_NODES=$(cat $BAD_NODES_FILE)
# echo "Note: bad nodes ($BAD_NODES) will be excluded"

# Slurm parameters
TASK=copy       # Task description
MEMO=4G          # Memory required
TIME=08:00:00 # Time required
CORE=1           # Cores required

# Assemble order prefix
ORDP="qsub -l h_rt="$TIME" -l h_vmem="$MEMO" -pe smp "$CORE" -R y -binding linear:"$CORE"" 

# Create directory for log files
LOGS=~/ukbiobank_knockoffs/data/logs
mkdir -p $LOGS

# Script to be run
SCRIPT="download_v2.sh"

# Define job name for this chromosome
JOBN=$TASK
OUTF=$LOGS"/"$JOBN".out"
ERRF=$LOGS"/"$JOBN".err"
# Assemble slurm order for this job
ORD=$ORDP" -N "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
# Print order
echo $ORD
# Submit order
$ORD
