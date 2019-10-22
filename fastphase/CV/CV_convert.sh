#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Convert data for CV.
#
# Author: Matteo Sesia
# Date:   07/17/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -n N -v V [-r] "
    echo "  -c chromosome   specify chromosome"
    echo "  -n samples      number of training samples"
    echo "  -v variants     number of variants in training set"
    echo "  -r              delete all intermediate results and recompute"
    exit 1
}

# Default values for optional input arguments
RESET=0
NTEST=1000

# If less than 1 arguments supplied, display usage
if [  $# -le 3 ]
then
  display_usage
  exit 1
fi

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":c:n:v:r" opt; do
  case $opt in
    c)
      echo "  - chromosome              : $OPTARG" >&2
      CHR=$OPTARG
      ;;
    n)
      echo "  - number of subjects      : $OPTARG" >&2
      NTRAIN=$OPTARG
      ;;
    v)
      echo "  - number of variants      : $OPTARG" >&2
      NVARIANTS=$OPTARG
      ;;
    r)
      echo "  - reset                   : TRUE" >&2
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
# Initialize variables
############################################################

# Load plink and datamash
export PATH=$PI_HOME/bin/:$PATH

# Initialize environment
DAT_DIR="/scratch/PI/candes/ukbiobank"               # Location of original data
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"           # Location of intermediate files

# Location of output
mkdir -p $TMP_DIR"/fastphase"
FP_DIR=$TMP_DIR"/fastphase/CV_"$NTRAIN
mkdir -p $FP_DIR
mkdir -p $FP_DIR"/input"

# Location of data for this chromosome
GENOTYPES=$DAT_DIR"/genotypes/ukb_gen_chr"$CHR
HAPLOTYPES=$DAT_DIR"/haplotypes/ukb_hap_chr"$CHR

# List of individuals
INDIVIDUALS=$TMP_DIR"/QC_output/individuals_QC.txt"

# Subsample individuals
INDIVIDUALS_SHUFFLE=$FP_DIR"/individuals.txt"
INDIVIDUALS_TRAIN=$FP_DIR"/individuals_train.txt"
INDIVIDUALS_TEST=$FP_DIR"/individuals_test.txt"
if [[ ! -f $INDIVIDUALS_TRAIN ]] || [[ ! -f $INDIVIDUALS_TEST ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Splitting observations into training/test sets"
  echo "----------------------------------------------------------------------------------------------------"
  sort -R $INDIVIDUALS > $INDIVIDUALS_SHUFFLE
  head -n $NTRAIN $INDIVIDUALS_SHUFFLE > $INDIVIDUALS_TRAIN  
  tail -n +$((NTRAIN+1)) $INDIVIDUALS_SHUFFLE | head -n $NTEST > $INDIVIDUALS_TEST
  rm $INDIVIDUALS_SHUFFLE
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping sample splitting because"
  echo $INDIVIDUALS_TRAIN" exists"
  echo $INDIVIDUALS_TEST" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
INDIVIDUALS=$INDIVIDUALS_SHUFFLE

# List of variants
VARIANTS=$TMP_DIR"/QC_output/QC_chr"$CHR".snplist"

# Subsample variants
if [ "$NVARIANTS" -gt "0" ]; then
  VARIANTS_SMALL=$FP_DIR"/input/QC_chr"$CHR".snplist"
  head -n $NVARIANTS $VARIANTS > $VARIANTS_SMALL
  VARIANTS=$VARIANTS_SMALL
fi

########################################
# Convert BGEN v1.2 into HAPS (train)
########################################
OUT_FILE=$FP_DIR"/input/ukb_hap_chr"$CHR"_train"
if [[ ! -f $OUT_FILE".haps" ]] || [[ ! -f $OUT_FILE".legend" ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting haplotypes on chromosome "$CHR" into HAPS format (train)"
  echo "----------------------------------------------------------------------------------------------------"
  plink2 --bgen $HAPLOTYPES".bgen" --sample $HAPLOTYPES".sample" --oxford-single-chr $CHR \
         --keep $INDIVIDUALS_TRAIN --extract $VARIANTS --export hapslegend --out $OUT_FILE
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes on chromosome "$CHR" into HAPS format because (train)"
  echo $OUT_FILE".haps exists"
  echo $OUT_FILE".legend exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

########################################
# Convert HAPS into INP (train)
########################################
OUT_FILE=$FP_DIR"/input/ukb_hap_chr"$CHR"_train"
if [[ ! -f $OUT_FILE".inp" ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting haplotypes on chromosome "$CHR" into INP format (train)"
  echo "----------------------------------------------------------------------------------------------------"
  ./haptoinp.sh $OUT_FILE
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes on chromosome "$CHR" into INP format (train) because"
  echo $OUT_FILE".inp exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

########################################
# Convert BGEN v1.2 into HAPS (test)
########################################
OUT_FILE=$FP_DIR"/input/ukb_hap_chr"$CHR"_test"
if [[ ! -f $OUT_FILE".haps" ]] || [[ ! -f $OUT_FILE".legend" ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting haplotypes on chromosome "$CHR" into HAPS format (test)"
  echo "----------------------------------------------------------------------------------------------------"
  plink2 --bgen $HAPLOTYPES".bgen" --sample $HAPLOTYPES".sample" --oxford-single-chr $CHR \
         --keep $INDIVIDUALS_TEST --extract $VARIANTS --export hapslegend --out $OUT_FILE
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes on chromosome "$CHR" into HAPS format (test) because"
  echo $OUT_FILE".haps exists"
  echo $OUT_FILE".legend exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

########################################
# Convert HAPS into INP (test)
########################################
OUT_FILE=$FP_DIR"/input/ukb_hap_chr"$CHR"_test"
if [[ ! -f $OUT_FILE".inp" ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Converting haplotypes on chromosome "$CHR" into INP format (test)"
  echo "----------------------------------------------------------------------------------------------------"
  ./haptoinp.sh $OUT_FILE
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping conversion of haplotypes on chromosome "$CHR" into INP format (test) because"
  echo $OUT_FILE".inp exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
