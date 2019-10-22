#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Impute missing values in test INP file
#
# Authors: Matteo Sesia
# Date:    07/20/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -k K [-u]"
    echo "  -c chromosome   specify chromosome"
    echo "  -k K            number of motifs in fastPHASE model"
    echo "  -u              whether to interpret the genotypes as unphased"
    exit 1
}

# Default values for optional input arguments
PHASE="phased"
PHASE_FLAG="-B"
RESET=0

# If less than 2 arguments supplied, display usage
if [  $# -le 3 ]
then
  display_usage
  exit 1
fi

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":c:k:u" opt; do
  case $opt in
    c)
      echo "  - chromosome            : $OPTARG" >&2
      CHR=$OPTARG
      ;;
    k)
      echo "  - fastPhase motifs      : $OPTARG" >&2
      K=$OPTARG
      ;;
    u)
      echo "  - unphased              : TRUE" >&2
      PHASE="unphased"
      PHASE_FLAG=""
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
# Impute masked values with fastPHASE
############################################################

# Load plink
export PATH=$PI_HOME/bin/:$PATH

# Location of input
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"
INP_DIR=$TMP_DIR"/fastphase/test"
INP_FILE=$INP_DIR/"ukb_hap_chr"$CHR"_masked.inp"
HMM_FILE=$TMP_DIR"/fastphase/"$PHASE"_K"$K"/ukb_hap_chr"$CHR

# Location of output
mkdir -p $TMP_DIR"/fastphase/test_imputed"
OUT_DIR=$TMP_DIR"/fastphase/test_imputed/"$PHASE"_K"$K
mkdir -p $OUT_DIR
OUT_FILE=$OUT_DIR"/ukb_hap_chr"$CHR

# Other fastPHASE parameters
NIT=0            # Number of iterations
SED=123          # Random seed

# Run fastPHASE on this chromosome
CMD="fastphase -Pp -T1 -K"$K" "$PHASE_FLAG" -C"$NIT" -S-7 -I"$HMM_FILE" -o"$OUT_FILE" "$INP_FILE
$CMD

# Convert imputed haplotypes into standard INP format
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

############################################################
# Evaluate prediction error
############################################################
echo "Evaluating imputation error..."
ml R/3.4.0
Rscript --vanilla test_eval.R $CHR $PHASE $K
