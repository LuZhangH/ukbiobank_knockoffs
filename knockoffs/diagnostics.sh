#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Compute knockoff diagnostics.
#
# Authors: Matteo Sesia
# Date:   12/13/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -m method -r factor -k K -s seed [-n subjects] [--reset]"
    echo "  -c chromosome         specify chromosome"
    echo "  -m clumping method    LD measure and clustering method (Radj, Rsingle, Dadj, Dsingle)"
    echo "  -r clumping factor    clumping percentage (0,100)"
    echo "  -k motifs             specify number of haoplotype motifs in the fastPHASE model"
    echo "  -s seed               random seed (default: 1)"
    echo "  --reset               delete all intermediate results and recompute"
    exit 1
}

# Set default options
RESET=0
METHOD="Rsingle"
SEED=1

# Call getopt to validate theinput.
options=$(getopt -o c:m:r:k:s: --long reset -- "$@")
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
    -m)
      shift;
      METHOD=$1
      [[ ! $METHOD =~ "Rsingle"|"Radj"|"Dsingle"|"Dadj" ]] && {
        echo "Incorrect value for option -m provided"
        exit 1
      }
      ;;
    -s)
      shift;
      SEED=$1
      ;;
    -r)
      shift;
      CLUMPING=$1
      ;;
    -k)
      shift;
      K=$1
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
if [[ -z "$CHR" ]] || [[ -z "$METHOD" ]] || [[ -z "$CLUMPING" ]] || [[ -z "$K" ]]; then
  display_usage
fi

# # Directories for input and output files
# TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"         # Location of intermediate files
# INP_DIR=$TMP_DIR"/knockoffs/"$METHOD$CLUMPING"_K"$K"_s"$SEED
# OUT_DIR=$TMP_DIR"/knockoff_diagnostics/"
# mkdir -p $OUT_DIR
# chmod go-rwx $OUT_DIR
# chmod g+rx $OUT_DIR
# OUT_DIR=$OUT_DIR$METHOD$CLUMPING"_K"$K"_s"$SEED
# mkdir -p $OUT_DIR
# chmod go-rwx $OUT_DIR
# chmod g+rx $OUT_DIR

# Directories for input and output files
DAT_DIR="/oak/stanford/groups/candes/ukbiobank_tmp"         # Location of intermediate files
INP_DIR=$DAT_DIR"/knockoffs/"$METHOD$CLUMPING"_K"$K
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"         # Location of intermediate files
OUT_DIR=$TMP_DIR"/knockoff_diagnostics/"
mkdir -p $OUT_DIR
chmod go-rwx $OUT_DIR
chmod g+rx $OUT_DIR
OUT_DIR=$OUT_DIR$METHOD$CLUMPING"_K"$K"_s"$SEED
mkdir -p $OUT_DIR
chmod go-rwx $OUT_DIR
chmod g+rx $OUT_DIR

# # Directories for input and output files
# TMP_DIR="/oak/stanford/groups/candes"
# INP_DIR=$TMP_DIR"/tmp/"$METHOD$CLUMPING"_K"$K
# OUT_DIR=$TMP_DIR"/tmp/"$METHOD$CLUMPING"_K"$K
# mkdir -p $OUT_DIR

# Location of this group's files
GENOTYPES=$INP_DIR"/ukb_gen_chr"$CHR
OUT_BASENAME=$OUT_DIR"/ukb_gen_chr"$CHR
FILE_1=$OUT_BASENAME".ld"
FILE_2=$OUT_BASENAME".frq"

if [[ ! -f $FILE_1 ]] || [[ ! -f $FILE_2 ]] || [[ $RESET -eq 1 ]]; then

  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo " Computing knockoff diagnostics for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo "----------------------------------------------------------------------------------------------------"
  plink --bfile $GENOTYPES \
        --freq \
        --r2 --ld-window 1000 --ld-window-kb 1000 --ld-window-r2 0.01 \
        --memory 9000 \
  	--out $OUT_BASENAME
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping knockoff diagnostics for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo $FILE_1" exists"
  echo $FILE_2" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
