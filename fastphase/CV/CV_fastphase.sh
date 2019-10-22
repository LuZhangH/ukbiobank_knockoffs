#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Run fastPHASE on masked data to tune K.
#
# Author: Matteo Sesia
# Date:   07/17/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -k K -n samples [-r] "
    echo "  -c chromosome   specify chromosome"
    echo "  -k K            specify fastPHASE model"
    echo "  -n N            number of individuals in training and test set"
    echo "  -r              delete all intermediate results and recompute"
    exit 1
}

# Default values for optional input arguments
RESET=0
PHASE="phased"
PHASE_FLAG="-B"

# If less than 1 arguments supplied, display usage
if [  $# -le 3 ]
then
  display_usage
  exit 1
fi

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":c:k:n:r" opt; do
  case $opt in
    c)
      echo "  - chromosome            : $OPTARG" >&2
      CHR=$OPTARG
      ;;
    k)
      echo "  - fastPhase motifs      : $OPTARG" >&2
      K=$OPTARG
      ;;
    n)
      echo "  - number of subjects      : $OPTARG" >&2
      NSUBJECTS=$OPTARG
      ;;
    r)
      echo "  - reset                 : TRUE" >&2
      RESET=1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

############################################################
# Run fastPHASE on training data
############################################################

# Load plink
export PATH=$PI_HOME/bin/:$PATH

# Location of input
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp/fastphase/CV_"$NSUBJECTS"/"
INP_FILE=$TMP_DIR"input/ukb_hap_chr"$CHR"_train.inp"

# Location of output
OUT_DIR=$TMP_DIR"K"$K
mkdir -p $OUT_DIR
OUT_FILE=$OUT_DIR"/ukb_hap_chr"$CHR

# Other fastPHASE parameters
NIT=20           # Number of iterations
SED=123          # Random seed
PHASE="phased"
PHASE_FLAG="-B"

if [[ ! -f $OUT_FILE"_finallikelihoods" ]] || [[ $RESET -eq 1 ]]; then
  # Run fastPHASE on this chromosome
  CMD="fastphase -Pp -T1 -g -H-4 -K"$K" "$PHASE_FLAG" -C"$NIT" -S"$SED" -o"$OUT_FILE" "$INP_FILE
  $CMD
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping model fitting for K=$K on chromosome $CHR because"
  echo $OUT_FILE"_finallikelihoods exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
OUT_FILE_TRAIN=$OUT_FILE

############################################################
# Run fastPHASE on test data
############################################################

# Location of input
INP_FILE=$TMP_DIR"input/mask01/ukb_hap_chr"$CHR"_test.inp"

# Location of output
OUT_DIR=$TMP_DIR"K"$K"/mask01"
mkdir -p $OUT_DIR
OUT_FILE=$OUT_DIR"/ukb_hap_chr"$CHR

# Other fastPHASE parameters
NIT=0           # Number of iterations
SED=123          # Random seed
PHASE="phased"
PHASE_FLAG="-B"

if [[ ! -f $OUT_FILE"_imputed.inp" ]] || [[ $RESET -eq 1 ]]; then
  # Run fastPHASE on this chromosome
  CMD="fastphase -Pp -T1 -K"$K" "$PHASE_FLAG" -C"$NIT" -S-7 -I"$OUT_FILE_TRAIN" -o"$OUT_FILE" "$INP_FILE
  $CMD
  # Convert output haplotype file into standard format
  # Remove blank lines
  HAP_FILE=$OUT_FILE"_haplotypes.out"
  awk 'NF' $HAP_FILE > $HAP_FILE".tmp"
  mv $HAP_FILE".tmp" $HAP_FILE
  # Add header
  HEADER=$(sed -n '2p' $INP_FILE)
  sed -i "1i "$HEADER $HAP_FILE
  HEADER=$(sed -n '1p' $INP_FILE)
  sed -i "1i "$HEADER $HAP_FILE
  # Rename file
  mv $HAP_FILE $OUT_FILE"_imputed.inp"
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping model fitting for K=$K on chromosome $CHR because"
  echo $OUT_FILE"_imputed.inp exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
