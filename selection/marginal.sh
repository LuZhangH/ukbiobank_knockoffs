#!/bin/bash
# UK Biobank GWAS
#
# Class: script
# 
# Compute marginal feature importance measures on variables and knockoffs.
#
# Authors: Matteo Sesia
# Date:  Sept 25, 2018

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -c chromosome -m clumping -k K -p phenotype [--reset --resample]"
    echo "  -c chromosome         specify chromosome"
    echo "  -m clumping method    LD measure, clustering method and clumping (Radj10, Radj50)"
    echo "  -k motifs             specify number of haoplotype motifs in the fastPHASE model"
    echo "  -p phenotype          specify number of haoplotype motifs in the fastPHASE model"
    echo "  -r regression         either linear or logistic"
    echo "  --reset               delete all intermediate results and recompute"
    exit 1
}

# Set default options
RESET=0
CLUMPING="Radj50"

# Call getopt to validate theinput.
options=$(getopt -o c:m:k:p:r: --long reset -- "$@")
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
      CLUMPING=$1
      ;;
    -k)
      shift;
      K=$1
      ;;
    -p)
      shift;
      PHENOTYPE=$1
      ;;
    -r)
      shift;
      REGRESSION=$1
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
if [[ -z "$CHR" ]] || [[ -z "$CLUMPING" ]] || [[ -z "$K" ]]; then
  display_usage
fi

# Load plink
export PATH=/share/PI/candes/bin/:$PATH

# Input genotypes
GENO_DIR="/scratch/PI/candes/ukbiobank_tmp/augmented_data/"$METHOD$CLUMPING"_K"$K
GENO_FILE=$GENO_DIR/"ukb_gen_chr"$CHR

# Input phenotypes
PHENO_FILE="/scratch/PI/candes/ukbiobank_tmp/association/phenotypes/phenotypes.tab"

###############
# Association #
###############

# Output file
OUT_DIR="/scratch/PI/candes/ukbiobank_tmp/association"
OUT_DIR=$OUT_DIR/$PHENOTYPE
mkdir -p $OUT_DIR
OUT_DIR=$OUT_DIR"/"$METHOD$CLUMPING"_K"$K
mkdir -p $OUT_DIR
OUT_FILE=$OUT_DIR"/ukb_chr"$CHR

# # Run association analysis
# plink --bfile $GENO_FILE \
#       --$REGRESSION hide-covar \
#       --pheno $PHENO_FILE \
#       --pheno-name $PHENOTYPE \
#       --covar $PHENO_FILE \
#       --covar-name PC.1 PC.2 PC.3 PC.4 PC.5 \
#       --out $OUT_FILE

# Rename output file
#mv $OUT_FILE".assoc."$REGRESSION $OUT_FILE".assoc"

############
# Clumping #
############
echo $OUT_FILE".assoc"
plink --bfile $GENO_FILE \
      --clump $OUT_FILE".assoc" \
      --clump-p1 5e-8 --clump-p2 1e-5 --clump-r2 0.1 --clump-kb 500 \
      --out $OUT_FILE
