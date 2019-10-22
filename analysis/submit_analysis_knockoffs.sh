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
#PHENOTYPE_LIST=("cvd" "diabetes" "hypothyroidism" "respiratory")
#PHENOTYPE_LIST=("glaucoma")
RESOLUTION_LIST=("Radj2" "Radj5" "Radj10" "Radj20" "Radj50" "Radj75" "Radj100")
#RESOLUTION_LIST=("Radj2")
SEED_LIST=("0")

# Slurm parameters
PART=candes,hns,owners               # Partition names
MEMO=5G                   # Memory required (20GB)
TIME=00-00:20:00          # Time required (12h)
CORE=1                    # Cores required (10)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME" --exclude=$BAD_NODES"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

# Loop over configurations and chromosomes
for SEED in "${SEED_LIST[@]}"; do
  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    for PHENOTYPE in "${PHENOTYPE_LIST[@]}"; do
      OUT_DIR="/scratch/PI/candes/ukbiobank_tmp/analysis/knockoffs"
      #OUT_FILE=$OUT_DIR"/"$PHENOTYPE"_"$RESOLUTION"_lasso_nopc.txt"
      OUT_FILE=$OUT_DIR"/"$PHENOTYPE"_"$RESOLUTION"_s"$SEED"_lasso_stats.txt"
      #if [ ! -f $OUT_FILE ]; then
        # Script to be run
        SCRIPT="analysis_knockoffs.sh $PHENOTYPE $RESOLUTION $SEED"
        # Define job name for this chromosome
        JOBN="knockoffs_"$PHENOTYPE"_"$RESOLUTION"_s"$SEED
        OUTF=$LOGS"/"$JOBN".out"
        ERRF=$LOGS"/"$JOBN".err"
        # Assemble slurm order for this job
        ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
        # Print order
        echo $ORD
        # Submit order
        $ORD
        # Run command now
        #./$SCRIPT
        #else 
        #echo $OUT_FILE" was found"
      #fi
    done
  done
done
