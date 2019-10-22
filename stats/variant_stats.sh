#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Convert phased haplotype files from BGEN 1.2 to VCF
# and then to the text file format required by fastPHASE. 
#
# Authors: Matteo Sesia
# Date:    08/08/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome [--reset]"
    echo "  -c chromosome   specify chromosome"
    echo "  --reset         delete all intermediate results and recompute"
    exit 1
}

# Set default options
RESET=0

# Call getopt to validate theinput.
options=$(getopt -o c: --long reset -- "$@")
[ $? -eq 0 ] || {
    echo "Incorrect options provided"
    exit 1
}
eval set -- "$options"
while true; do
    case "$1" in
    -c)
      shift;
      CHR=$1
      [[ ! $CHR =~ 1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22 ]] && {
        echo "Incorrect value for option -c provided"
        exit 1
      }
      ;;
    --reset)
        RESET=1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done

# Verify that all required parameters have been supplied
if [[ -z "$CHR" ]]; then
  display_usage
fi


############################################################
# Initialize variables
############################################################

# Load plink and datamash
export PATH=$PI_HOME/bin/:$PATH

# Initialize environment
DAT_DIR="/scratch/PI/candes/ukbiobank"               # Location of original data
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"           # Location of intermediate files

# Location of output
OUT_DIR=$TMP_DIR"/stats"
mkdir -p $OUT_DIR
OUT_BASENAME=$OUT_DIR"/ukb_gen_chr"$CHR
chmod ugo-rwx $OUT_DIR
chmod u+rwx $OUT_DIR

# Location of data for this chromosome
GENOTYPES=$DAT_DIR"/genotypes/ukb_gen_chr"$CHR
HAPLOTYPES=$DAT_DIR"/haplotypes/ukb_hap_chr"$CHR

# List of individuals and variants that passed QC
INDIVIDUALS=$TMP_DIR"/QC_output/individuals_QC.txt"
VARIANTS=$TMP_DIR"/QC_output/QC_chr"$CHR".snplist"

##############################
# Frequency
##############################
OUT_FILE=$OUT_BASENAME".frq"
# Compute variant frequencies if not already present
if [[ ! -f $OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing variant frequencies for chromosome "$CHR
  echo "----------------------------------------------------------------------------------------------------"
  plink --bfile $GENOTYPES --keep $INDIVIDUALS --extract $VARIANTS \
        --freq \
	--out $OUT_BASENAME
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping variant frequencies for chromosome "$CHR" because"
  echo $OUT_FILE" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

##############################
# Correlation
##############################
OUT_BASENAME=$OUT_DIR"/ukb_gen_chr"$CHR
OUT_FILE=$OUT_BASENAME".ld"
# Compute correlations if not already present
if [[ ! -f $OUT_FILE ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing correlation matrix for chromosome "$CHR
  echo "----------------------------------------------------------------------------------------------------"
  plink --bfile $GENOTYPES --keep $INDIVIDUALS --extract $VARIANTS \
        --r2 dprime --ld-window 1000 --ld-window-kb 1000 --ld-window-r2 0.01 \
	--out $OUT_BASENAME
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping correlation matrix for chromosome "$CHR" because"
  echo $OUT_FILE" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
