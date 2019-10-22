#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Export augmented genetic data to shared folder.
#
# Authors: Matteo Sesia
# Date:   08/09/2018

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
    -r)
      shift;
      CLUMPING=$1
      ;;
    -k)
      shift;
      K=$1
      ;;
    -s)
      shift;
      SEED=$1
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

# Initialize environment
DAT_DIR="/scratch/PI/candes/ukbiobank"             # Location of original data
TMP_DIR="/scratch/PI/candes/ukbiobank_tmp"         # Location of intermediate files

##############################
# Copy the encrypted data to a shared folder
##############################
# Directory for output files
INP_DIR=$TMP_DIR"/knockoffs/"$METHOD$CLUMPING"_K"$K"_s"$SEED
OUT_DIR=$TMP_DIR"/augmented_data/"
mkdir -p $OUT_DIR
chmod ugo-rwx $OUT_DIR
chmod g+rx $OUT_DIR
chmod u+rwx $OUT_DIR
OUT_DIR=$OUT_DIR$METHOD$CLUMPING"_K"$K"_s"$SEED
mkdir -p $OUT_DIR
chmod ugo-rwx $OUT_DIR
chmod g+rx $OUT_DIR
chmod u+rwx $OUT_DIR

# Location of this group's files
FILE_1="/ukb_gen_chr"$CHR".bed"
FILE_2="/ukb_gen_chr"$CHR".bim"
FILE_3="/ukb_gen_chr"$CHR".fam"
FILE_4_IN=$TMP_DIR"/clumping/"$METHOD$CLUMPING"/grp_chr"$CHR".txt"
FILE_4="/ukb_gen_chr"$CHR".grp"

if [[ ! -f $OUT_DIR$FILE_1 ]] || [[ ! -f $OUT_DIR$FILE_2 ]] || [[ ! -f $OUT_DIR$FILE_3 ]] || [[ ! -f $OUT_DIR$FILE_4 ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo " Exporting augmented genotype data for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo "----------------------------------------------------------------------------------------------------"
  rsync -auv $INP_DIR$FILE_1 $OUT_DIR$FILE_1
  rsync -auv $INP_DIR$FILE_2 $OUT_DIR$FILE_2
  rsync -auv $INP_DIR$FILE_3 $OUT_DIR$FILE_3
  rsync -auv $FILE_4_IN $OUT_DIR$FILE_4
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping export of augmented genotype data for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo $OUT_DIR$FILE_1" exists"
  echo $OUT_DIR$FILE_2" exists"
  echo $OUT_DIR$FILE_3" exists"
  echo $OUT_DIR$FILE_4" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi
