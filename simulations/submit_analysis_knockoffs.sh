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
RESOLUTION_LIST=("Radj2" "Radj5" "Radj10" "Radj20" "Radj50" "Radj75" "Radj100")
#RESOLUTION_LIST=("Radj75" "Radj100")
FOLD_LIST=$(seq -f "%g" 1 1)

# Slurm parameters
PART=candes,hns,owners       # Partition names
#MEMO=20G                 # Memory required
#TIME=01-00:00:00          # Time required (8h)
#CORE=10                   # Cores required
MEMO=5G
TIME=00-00:20:00
CORE=1

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Loop over configurations and chromosomes
for FOLD in $FOLD_LIST; do
  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    cat parameters_p1.txt | while read PARAM_LINE
    do 
      PARAM_LINE=($PARAM_LINE)
      EXPERIMENT_NAME=${PARAM_LINE[0]}
      EXPERIMENT_NUMB=${PARAM_LINE[1]}

      OUT_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations/lasso/fold_"$FOLD"/"$RESOLUTION
      OUT_FILE=$OUT_DIR"/lasso_"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB".txt2"
      
      if [ ! -f $OUT_FILE ]; then
        echo "Preparing experiment: ("$EXPERIMENT_NAME","$EXPERIMENT_NUMB"), resolution: "$RESOLUTION", fold: "$FOLD
        # Script to be run
        SCRIPT="analysis_knockoffs.sh $EXPERIMENT_NAME $EXPERIMENT_NUMB $FOLD $RESOLUTION"
        # Define job name for this chromosome
        JOBN="knockoffs_"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB"_"$FOLD"_"$RESOLUTION
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
done
