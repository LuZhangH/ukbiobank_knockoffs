#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Perform clumping, generate knockoffs, augment the genetic data and export it to the shared folder.
#
# Authors: Matteo Sesia
# Date:    08/09/2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -m method -r factor -k K [-n subjects] [--reset --resample]"
    echo "  -c chromosome         specify chromosome"
    echo "  -m clumping method    LD measure and clustering method (Radj, Rsingle, Dadj, Dsingle)"
    echo "  -r clumping factor    clumping percentage (0,100)"
    echo "  -k motifs             specify number of haoplotype motifs in the fastPHASE model"
    echo "  -n subjects           number of subjects to construct knockoffs for (default: all)"
    echo "  -s seed               random seed (default: 1)"
    echo "  --gaussian            solve the SDP for Gaussian knockoffs"
    echo "  --reset               delete all intermediate results and recompute"
    echo "  --resample            delete knockoffs and resample"
    exit 1
}

# Set default options
RESET=0
RESAMPLE=0
METHOD="Rsingle"
NSUBJECTS=-1
GAUSSIAN=0
SEED=1

# Call getopt to validate theinput.
options=$(getopt -o c:m:r:k:n:s: --long reset,resample,gaussian -- "$@")
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
    -n)
      shift;
      NSUBJECTS=$1
      ;;
    -s)
      shift;
      SEED=$1
      ;;
    --reset)
        RESET=1
        ;;
    --gaussian)
        GAUSSIAN=1
        ;;
    --resample)
        RESAMPLE=1
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
# Clumping
##############################

ml gcc
ml R/3.5.1
export PATH=$PI_HOME/bin/:$PATH                    # Load PLINK

# Location of data for this chromosome
GENOTYPES=$DAT_DIR"/genotypes/ukb_gen_chr"$CHR
HAPLOTYPES=$DAT_DIR"/haplotypes/ukb_hap_chr"$CHR

# List of individuals and variants that were fed in to fastPHASE
INDIVIDUALS=$TMP_DIR"/fastphase/data/individuals.txt"
VARIANTS=$TMP_DIR"/fastphase/data/ukb_hap_chr"$CHR".legend"

# Create directory for clumping files
GRP_DIR=$TMP_DIR"/clumping"
mkdir -p $GRP_DIR

# Create directory for dendrograms
mkdir -p $GRP_DIR"/"$METHOD

# Create directory for group files
OUT_DIR=$GRP_DIR"/"$METHOD$CLUMPING
mkdir -p $OUT_DIR

# Location of group files
OUT_FILE_1=$OUT_DIR"/grp_chr"$CHR".txt"
OUT_FILE_2=$OUT_DIR"/rep_chr"$CHR".txt"

# Compute variant frequencies if not already present
if [[ ! -f $OUT_FILE_1 ]] || [[ ! -f $OUT_FILE_2 ]] || [[ $RESET -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Making groups for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo "----------------------------------------------------------------------------------------------------"
  Rscript --vanilla clumping.R $CHR $CLUMPING $METHOD
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping groups for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD" because"
  echo $OUT_FILE_1" exists"
  echo $OUT_FILE_2" exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

##############################
# HMM knockoffs
##############################
# Create directory for knockoff files
KNO_DIR=$TMP_DIR"/knockoffs"
mkdir -p $KNO_DIR
OUT_DIR=$KNO_DIR"/"$METHOD$CLUMPING"_K"$K"_s"$SEED
mkdir -p $OUT_DIR

# Location of this group's files
OUT_FILE=$OUT_DIR"/ukb_gen_chr"$CHR

# Generate knockoffs if not already present
if [[ ! -f $OUT_FILE".bed" ]] || [[ ! -f $OUT_FILE".bim" ]] || [[ ! -f $OUT_FILE".fam" ]] || [[ ! -f $OUT_FILE".key" ]] || [[ $RESET -eq 1 ]] || [[ $RESAMPLE -eq 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Making knockoffs for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo "----------------------------------------------------------------------------------------------------"
  # Generate knockoffs and save them in PED format
  Rscript --vanilla knockoffs.R $CHR $CLUMPING $METHOD $K $NSUBJECTS $SEED
  # Convert PED to BED
  plink --file $OUT_FILE --no-fid --no-parents --no-sex --no-pheno --keep-allele-order --make-bed --out $OUT_FILE
  # Add subject information to FAM file  
  SEX_FILE=$TMP_DIR"/fastphase/data/ukb_hap_chr"$CHR".sample"  
  awk 'FNR > 2 { print $1,$4 }' $SEX_FILE > $OUT_FILE".sex"
  awk '{ print $1,$2,$3,$4,$6 }' $OUT_FILE".fam"> $OUT_FILE".fam.tmp"
  join -1 1 -2 1 $OUT_FILE".fam.tmp" $OUT_FILE".sex" > $OUT_FILE".fam.tmp2"
  awk '{ print $1,$2,$3,$4,$6,$5 }' $OUT_FILE".fam.tmp2"> $OUT_FILE".fam"
  rm $OUT_FILE".fam.tmp"
  rm $OUT_FILE".fam.tmp2"
  rm $OUT_FILE".sex"
  # Verify that the knockoffs were constructed correctly
  #Rscript --vanilla verify.R $CHR $CLUMPING $METHOD $K 2000
  # Remove PED file
  rm $OUT_FILE".ped"
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping knockoffs for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo $OUT_FILE".bed exists"
  echo $OUT_FILE".bim exists"
  echo $OUT_FILE".fam exists"
  echo $OUT_FILE".key exists"
  echo "----------------------------------------------------------------------------------------------------"
fi

##############################
# Gaussian knockoffs
##############################
# Create directory for knockoff files
KNO_DIR=$TMP_DIR"/knockoffs_gau"
mkdir -p $KNO_DIR

# Create directory for knockoff files
OUT_DIR=$KNO_DIR"/"$METHOD$CLUMPING
mkdir -p $OUT_DIR

# Location of this group's files
OUT_FILE=$OUT_DIR"/ukb_gen_chr"$CHR".diags"

# Generate Gaussian knockoffs
if [[ $GAUSSIAN -eq 1 ]] && ([[ ! -f $OUT_FILE ]] || [[ $RESET -eq 1 ]]) ; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Making gaussian knockoffs for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo "----------------------------------------------------------------------------------------------------"
  Rscript --vanilla knockoffs_gaussian.R $CHR $CLUMPING $METHOD $K
else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping gaussian knockoffs for chromosome "$CHR" with "$CLUMPING"% clumping and method "$METHOD
  echo $OUT_FILE" exists or --gaussian option was not set"
  echo "----------------------------------------------------------------------------------------------------"
fi
