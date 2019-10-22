#!/bin/bash

################
# Parse input  #
################

# EXPERIMENT_NAME="six"
# EXPERIMENT_NUMB="10"
# FOLD="1"

EXPERIMENT_NAME=$1
EXPERIMENT_NUMB=$2
FOLD=$3
RESOLUTION=$4

#########
# Setup #
#########

# Whether to run LASSO
RUN_LASSO=1

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

echo "----------------------------------------------------------------------------------------------------"
echo "Performing association analysis with knockoffs using the following parameters:"
echo "  - Experiment        :" $EXPERIMENT_NAME
echo "  - Experiment number :" $EXPERIMENT_NUMB
echo "  - Fold              :" $FOLD
echo "  - Resolution        :" $RESOLUTION
echo "  - Chromosome range  :" $CHR_MIN".."$CHR_MAX
echo "  - Run-LASSO         :" $RUN_LASSO
echo "----------------------------------------------------------------------------------------------------"
echo ""

#################################
# Initialize input/output paths #
#################################

# Input
PHENO_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small/phenotypes"

# Output
STATS_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small/bolt_screen"
STATS_FILE=$STATS_DIR"/fold_"$FOLD"/stats_"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB".txt"

if [ $RUN_LASSO -eq 1 ]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing importance measures with LASSO"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
  LASSO_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations_small/lasso"
  mkdir -p $LASSO_DIR
  mkdir -p $LASSO_DIR"/fold_"$FOLD
  mkdir -p $LASSO_DIR"/fold_"$FOLD"/"$RESOLUTION
  source activate ukb
  Rscript --vanilla lasso.R $EXPERIMENT_NAME $EXPERIMENT_NUMB $FOLD $RESOLUTION
else
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping computation of importance measures with LASSO"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
fi
