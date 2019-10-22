#!/bin/bash

################
# Parse input  #
################

PHENOTYPE=$1

#########
# Setup #
#########

# Whether to fit BOLT-LMM
RUN_BOLT=1
# Whether to perform clumping with PLINK
RUN_PLINK=1

# List of thresholds for clumping
CLUMP_THRESHOLD_LIST=(0.00000005)

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)
FOLD="01"

# Number of threads
THREADS=10

echo "----------------------------------------------------------------------------------------------------"
echo "Performing association analysis with BOLT-LMM using the following parameters:"
echo "  - Phenotype         :" $PHENOTYPE
echo "  - Chromosome range  :" $CHR_MIN".."$CHR_MAX
echo "  - FOLD              :" $FOLD
echo "  - Run-BOLT          :" $RUN_BOLT
echo "  - Run-PLINK         :" $RUN_PLINK
echo "  - Run-CAVIAR        :" $RUN_CAVIAR
echo "----------------------------------------------------------------------------------------------------"
echo ""

#################################
# Initialize input/output paths #
#################################

# Stuff for bolt
LD_TABLE="/home/groups/candes/Software/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz"
MAP_TABLE="/home/groups/candes/Software/BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz"

# Input
PHENO_DIR="/scratch/groups/candes/ukbiobank_tmp/phenotypes"
PHENO_FILE=$PHENO_DIR"/phenotypes_qc.tab"
GENO_FILE="/scratch/groups/candes/ukbiobank/genotypes/ukb_gen_chr"

# List of excluded SNPs
SNP_DIR="/scratch/groups/candes/ukbiobank_tmp/analysis/variants"
SNP_EXCLUDE_FILE=$SNP_DIR"/ukb_gen_exclude.txt"

# List of subjects
FAM_FILE="/scratch/groups/candes/ukbiobank_tmp/analysis/individuals/ukb_individuals.fam"
FOLD_FILE="/scratch/PI/candes/ukbiobank_tmp/analysis/individuals_split/ukb_gen_"$FOLD".fam"
FOLD_EXCLUDE_FILE="/scratch/PI/candes/ukbiobank_tmp/analysis/individuals_split/ukb_gen_"$FOLD"_exclude.fam"

# Output
mkdir -p "/scratch/groups/candes/ukbiobank_tmp/analysis"
mkdir -p "/scratch/groups/candes/ukbiobank_tmp/analysis/"
STATS_DIR="/scratch/groups/candes/ukbiobank_tmp/analysis/bolt"
mkdir -p $STATS_DIR
STATS_FILE=$STATS_DIR"/"$PHENOTYPE"_stats_fold_"$FOLD".txt"

############
# BOLT-LLM #
############

if [ $RUN_BOLT -eq 1 ]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing association statistics with BOLT-LMM"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
 
  # Read list of covariates for BOLT
  COVARIATES_FILE=$STATS_DIR"/"$PHENOTYPE"_covariates.txt"
  COVARIATE_LIST=$(<$COVARIATES_FILE)
  echo "Covariates for "$PHENOTYPE":"
  echo $COVARIATE_LIST

  CALL_BOLT="
  bolt \
    --bed=$GENO_FILE"{$CHR_MIN:$CHR_MAX}.bed" \
    --bim=$GENO_FILE"{$CHR_MIN:$CHR_MAX}.bim" \
    --fam=$FAM_FILE \
    --remove=$FOLD_EXCLUDE_FILE \
    --exclude=$SNP_EXCLUDE_FILE \
    --maxMissingPerSnp=1 \
    --phenoFile=$PHENO_FILE \
    --phenoCol=$PHENOTYPE \
    --covarFile=$PHENO_FILE \
    $COVARIATE_LIST    
    --LDscoresFile=$LD_TABLE \
    --LDscoresMatchBp \
    --geneticMapFile=$MAP_TABLE \
    --statsFile=$STATS_FILE \
    --lmm \
    --numThreads=$THREADS
  "
  echo "Running BOLT with the following command:"
  echo $CALL_BOLT
  $CALL_BOLT
  echo "Output file: "$STATS_FILE
else
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping association statistics with BOLT-LMM"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
fi

##################################
# Clumping and marginal p-values #
##################################

if [ $RUN_PLINK -eq 1 ]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Clumping LMM discoveries with PLINK"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""

  # Output directory for clumping
  CLUMP_DIR="/scratch/groups/candes/ukbiobank_tmp/analysis/bolt_clumped"
  mkdir -p $CLUMP_DIR

  # Create temporary stats file with modified header
  STATS_FILE_TMP=$CLUMP_DIR"/"$PHENOTYPE"_stats_fold_"$FOLD".tmp"
  head -n 1 $STATS_FILE > $STATS_FILE_TMP
  sed -i 's/\S\+$/P/' $STATS_FILE_TMP
  tail -n +2 $STATS_FILE >> $STATS_FILE_TMP

  for CHR in $CHR_LIST; do
    # Clumping of LMM p-values
    for CLUMP_THRESHOLD in "${CLUMP_THRESHOLD_LIST[@]}"; do
      CLUMP_FILE=$CLUMP_DIR"/"$PHENOTYPE"_fold_"$FOLD"_"$CLUMP_THRESHOLD
      if [ ! -f $CLUMP_FILE"_chr"$CHR".clumped" ]; then
        echo "Clumping chromosome $CHR at threshold $CLUMP_THRESHOLD"
        plink --bfile $GENO_FILE$CHR \
              --exclude $SNP_EXCLUDE_FILE \
              --clump $STATS_FILE_TMP \
              --clump-p1 $CLUMP_THRESHOLD \
              --clump-r2 0.01 \
              --clump-kb 5000 \
              --out $CLUMP_FILE"_chr"$CHR
        rm $CLUMP_FILE"_chr"$CHR".log"
        rm $CLUMP_FILE"_chr"$CHR".nosex"
      else
        echo "Skipping clumping of chromosome $CHR at threshold "$CLUMP_THRESHOLD
      fi
    done
  done

  #Clumping after BHq
  echo "Applying BH for FDR at level 0.1, then clumping"
  Rscript --vanilla lmm_bh.R $PHENOTYPE

  # Erase temporary stats file with modified header
  rm $STATS_FILE_TMP

  CLUMP_THRESHOLD_LIST_BH=(BH ${CLUMP_THRESHOLD_LIST[@]})
  #Combine clumped files
  for CLUMP_THRESHOLD in "${CLUMP_THRESHOLD_LIST_BH[@]}"; do
    echo "Combining clumping files at threshold $CLUMP_THRESHOLD"
    CLUMP_FILE=$CLUMP_DIR"/"$PHENOTYPE"_fold_"$FOLD"_"$CLUMP_THRESHOLD
    # Write header
    HEADER=" CHR    F             SNP         BP        P    TOTAL   NSIG    S05    S01   S001  S0001    SP2"
    echo $HEADER > $CLUMP_FILE".tab"
    # Combine results of clumping into a single file
    for CHR in $(seq $CHR_MIN $CHR_MAX); do
      ls $CLUMP_FILE"_chr"$CHR".clumped"
      if [ -s $CLUMP_FILE"_chr"$CHR".clumped" ]; then
        tail -n +2 $CLUMP_FILE"_chr"$CHR".clumped" >> $CLUMP_FILE".tab"
        rm $CLUMP_FILE"_chr"$CHR".clumped"
      fi
    done
    # Remove empty lines
    sed -i '/^$/d' $CLUMP_FILE".tab"
    echo "Results written on "$CLUMP_FILE".tab"
  done

  # Parse clumped p-values and summarise discoveries  
  source activate ukb
  for CLUMP_THRESHOLD in "${CLUMP_THRESHOLD_LIST_BH[@]}"; do
    echo "Summarizing discoveries with clumping "$CLUMP_THRESHOLD
    Rscript --vanilla summarise_lmm.R $PHENOTYPE $CLUMP_THRESHOLD $FOLD
    SUMMARY_DIR="/scratch/groups/candes/ukbiobank_tmp/discoveries"
    SUMMARY_VARIANTS=$SUMMARY_DIR"/"$PHENOTYPE"_lmm_variants_fold_"$FOLD".txt"
    SUMMARY_REGIONS=$SUMMARY_DIR"/"$PHENOTYPE"_lmm_regions_fold_"$FOLD".txt"
    if [ $CLUMP_THRESHOLD == "BH" ]; then
      SUMMARY_VARIANTS_BH=$SUMMARY_DIR"/"$PHENOTYPE"_lmm_variants_fold_"$FOLD"_BH.txt"
      SUMMARY_REGIONS_BH=$SUMMARY_DIR"/"$PHENOTYPE"_lmm_regions_fold_"$FOLD"_BH.txt"
      mv $SUMMARY_VARIANTS $SUMMARY_VARIANTS_BH
      mv $SUMMARY_REGIONS $SUMMARY_REGIONS_BH
      SUMMARY_VARIANTS=$SUMMARY_VARIANTS_BH
      SUMMARY_REGIONS=$SUMMARY_VARIANTS_BH
    fi
    echo "Output written to:"
    ls $SUMMARY_VARIANTS
    ls $SUMMARY_REGIONS
  done

else
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping marginal p-values and clumping of LMM statistics with PLINK"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""

fi
