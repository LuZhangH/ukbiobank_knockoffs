#!/bin/bash

################
# Parse input  #
################

# EXPERIMENT_NAME="five"
# EXPERIMENT_NUMB="10"
# FOLD="1"

EXPERIMENT_NAME=$1
EXPERIMENT_NUMB=$2
FOLD=$3

#########
# Setup #
#########

# Whether to prune variables before fitting BOLT-LMM
RUN_PRUNING=1
# Whether to fit BOLT-LMM (as opposed to clumping only)
RUN_BOLT=1
# Whether to perform clumping with PLINK
RUN_PLINK=0
# Whether to perform fine mapping with CAVIAR
RUN_CAVIAR=0

# List of thresholds for clumping
CLUMP_THRESHOLD_LIST=(0.00000005 0.0000005 0.000005 0.00001 0.00005)

# Range of chromosomes to include in the analysis
CHR_MIN=1
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

# Number of threads
THREADS=10

echo "----------------------------------------------------------------------------------------------------"
echo "Performing association analysis with BOLT-LMM using the following parameters:"
echo "  - Experiment        :" $EXPERIMENT_NAME
echo "  - Experiment number :" $EXPERIMENT_NUMB
echo "  - Fold              :" $FOLD
echo "  - Chromosome range  :" $CHR_MIN".."$CHR_MAX
echo "  - Run-PRUNING       :" $RUN_PRUNING
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
PHENO_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/phenotypes"
PHENO_FILE=$PHENO_DIR"/"$EXPERIMENT_NAME"_phenotypes.tab"
GENO_FILE="/scratch/groups/candes/ukbiobank/genotypes/ukb_gen_chr"
PHENOTYPE="Y_a"$EXPERIMENT_NUMB"_"$FOLD

# List of subjects
FAM_FILE="/scratch/groups/candes/ukbiobank_tmp/simulations/subjects/ukb_gen.fam"
REMOVE_FILE="/scratch/groups/candes/ukbiobank_tmp/simulations/subjects/ukb_gen_remove.txt"

# List of excluded SNPs
SNP_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/variants"
SNP_EXCLUDE_FILE=$SNP_DIR"/ukb_gen_exclude.txt"

# Screening for BOLT
BOLT_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/bolt"
mkdir -p $BOLT_DIR
mkdir -p $BOLT_DIR"/fold_"$FOLD
mkdir -p $BOLT_DIR"/fold_"$FOLD"/screened"
BOLT_MODEL=$BOLT_DIR"/fold_"$FOLD"/screened/model_"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB".txt"

# Output
STATS_DIR="/scratch/groups/candes/ukbiobank_tmp/simulations/bolt/fold_"$FOLD
STATS_FILE=$STATS_DIR"/stats_"$EXPERIMENT_NAME"_"$EXPERIMENT_NUMB".txt"

######################################
# Prune SNPs before calling BOLT-LMM #
######################################

pruneSNPs(){
  plink \
    --bed $GENO_FILE"$1.bed" \
    --bim $GENO_FILE"$1.bim" \
    --fam $GENO_FILE"$1.fam" \
    --remove $REMOVE_FILE \
    --exclude $SNP_EXCLUDE_FILE \
    --indep-pairwise 50 100 0.5 \
    --out $BOLT_MODEL"_chr"$1
  rm $BOLT_MODEL"_chr"$1".log"
  rm $BOLT_MODEL"_chr"$1".nosex"
}

if [ $RUN_PRUNING -eq 1 ] && [ ! -f $BOLT_MODEL".pruned" ]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Pruning variants with PLINK before calling BOLT-LMM"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
  # Prune list of filtered SNPs for each chromosome
  for CHR in $CHR_LIST; do
    ((i=i%THREADS)); ((i++==0)) && wait
    pruneSNPs "$CHR" &
  done
  wait

  # Combine lists of pruned SNPs
  if [ -f $BOLT_MODEL".pruned" ] ; then
    rm $BOLT_MODEL".pruned"
  fi
  touch $BOLT_MODEL".pruned"
  for CHR in $CHR_LIST; do
    cat $BOLT_MODEL"_chr"$CHR".prune.in" >> $BOLT_MODEL".pruned"
    rm $BOLT_MODEL"_chr"$CHR".prune.in" $BOLT_MODEL"_chr"$CHR".prune.out"
    rm $BOLT_MODEL"_chr"$CHR".log"
  done
else
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping pruning of variants with PLINK before calling BOLT-LMM"
  echo $BOLT_MODEL".pruned exists"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
fi
BOLT_MODEL=$BOLT_MODEL".pruned"

############
# BOLT-LLM #
############

if [ $RUN_BOLT -eq 1 ]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Computing association statistics with BOLT-LMM"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
  bolt \
    --bed=$GENO_FILE"{$CHR_MIN:$CHR_MAX}.bed" \
    --bim=$GENO_FILE"{$CHR_MIN:$CHR_MAX}.bim" \
    --fam=$FAM_FILE \
    --remove=$REMOVE_FILE \
    --exclude=$SNP_EXCLUDE_FILE \
    --modelSnps=$BOLT_MODEL \
    --maxMissingPerSnp=1 \
    --phenoFile=$PHENO_FILE \
    --phenoCol=$PHENOTYPE \
    --LDscoresFile=$LD_TABLE \
    --LDscoresMatchBp \
    --geneticMapFile=$MAP_TABLE \
    --statsFile=$STATS_FILE \
    --lmm \
    --numThreads=$THREADS
  echo "Output file: "$STATS_FILE

else
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping association statistics with BOLT-LMM"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
fi

exit

##################################
# Clumping and marginal p-values #
##################################

if [ $RUN_PLINK -eq 1 ]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Marginal p-values and clumping LMM statistics with PLINK"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""

  # Output directories for marginal p-values
  MARGINAL_DIR="/scratch/PI/candes/ukbiobank_tmp/simulations/pvalues"
  mkdir -p $MARGINAL_DIR
  MARGINAL_DIR=$MARGINAL_DIR"/"$EXPERIMENT
  mkdir -p $MARGINAL_DIR
  MARGINAL_DIR=$MARGINAL_DIR"/fold_"$FOLD
  mkdir -p $MARGINAL_DIR

  # Output directory for clumping
  CLUMP_DIR=$STATS_DIR"/clumped"
  mkdir -p $CLUMP_DIR

  # # Heal p-values file (this will need to be removed later)
  # head -n 1 $STATS_FILE > $STATS_FILE".tmp"  
  # sed -i "s/P	P_BOLT_LMM/P_INF	P_BOLT_LMM/g" $STATS_FILE".tmp"
  # tail -n +2 $STATS_FILE >> $STATS_FILE".tmp"
  # mv $STATS_FILE".tmp" $STATS_FILE

  # Create temporary stats file with modified header
  head -n 1 $STATS_FILE > $STATS_FILE".tmp"
  sed -i 's/\S\+$/P/' $STATS_FILE".tmp"
  tail -n +2 $STATS_FILE >> $STATS_FILE".tmp"

  for CHR in $CHR_LIST; do
    # Marginal p-values
    MARGINAL_FILE=$MARGINAL_DIR"/s"$N_SIGNALS"_"$PHENOTYPE"_chr"$CHR
    if [ ! -s $MARGINAL_FILE".assoc.linear" ]; then
      echo "Computing marginal p-values for chromosome $CHR"   
      plink --bfile $GENO_FILE$CHR \
            --exclude $SNP_EXCLUDE_FILE \
            --linear \
            --pheno $PHENO_FILE \
            --pheno-name $PHENOTYPE \
            --out $MARGINAL_FILE
      rm $MARGINAL_FILE".log"
      rm $MARGINAL_FILE".nosex"
    else
      echo "Skipping marginal p-values for chromosome $CHR"   
    fi

    # Clumping of LMM p-values
    for CLUMP_THRESHOLD in "${CLUMP_THRESHOLD_LIST[@]}"; do
      CLUMP_FILE=$CLUMP_DIR"/stats_"$EXPERIMENT"_s"$N_SIGNALS"_"$PHENOTYPE"_"$CLUMP_THRESHOLD
      if [ ! -f $CLUMP_FILE"_chr"$CHR".clumped" ]; then
        echo "Clumping chromosome $CHR at threshold $CLUMP_THRESHOLD"
        plink --bfile $GENO_FILE$CHR \
              --exclude $SNP_EXCLUDE_FILE \
              --clump $STATS_FILE".tmp" \
              --clump-p1 $CLUMP_THRESHOLD \
              --clump-r2 0.01 \
              --clump-kb 5000 \
              --out $CLUMP_FILE"_chr"$CHR
        rm $CLUMP_FILE"_chr"$CHR".log"
        rm $CLUMP_FILE"_chr"$CHR".nosex"
      else
        echo "Skipping clumping of chromosome $CHR at threshold $CLUMP_THRESHOLD"
      fi
    done
  done

  # Erase temporary stats file with modified header
  rm $STATS_FILE".tmp"

  for CLUMP_THRESHOLD in "${CLUMP_THRESHOLD_LIST[@]}"; do
    echo "Combining clumping files at threshold $CLUMP_THRESHOLD"
    CLUMP_FILE=$CLUMP_DIR"/stats_"$EXPERIMENT"_s"$N_SIGNALS"_"$PHENOTYPE"_"$CLUMP_THRESHOLD
    # Write header
    HEADER=" CHR    F             SNP         BP        P    TOTAL   NSIG    S05    S01   S001  S0001    SP2"
    echo $HEADER > $CLUMP_FILE".clumped"
    # Combine results of clumping into a single file
    for CHR in $(seq 1 22); do
      tail -n +2 $CLUMP_FILE"_chr"$CHR".clumped" >> $CLUMP_FILE".clumped"
    done
    echo "Results written on $CLUMP_FILE.clumped"
  done

else
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping marginal p-values and clumping of LMM statistics with PLINK"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""

fi

################
# Fine mapping #
################

if [ $RUN_CAVIAR -eq 1 ]; then
  echo "----------------------------------------------------------------------------------------------------"
  echo "Performing fine mapping with CAVIAR"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
  mkdir -p "/scratch/PI/candes/ukbiobank_tmp/simulations/caviar/"
  for CLUMP_THRESHOLD in "${CLUMP_THRESHOLD_LIST[@]}"; do
    echo "Applying fine mapping to results clumped at threshold "$CLUMP_THRESHOLD
    Rscript --vanilla fine_mapping.R $EXPERIMENT $N_SIGNALS $PHENOTYPE $FOLD $CLUMP_THRESHOLD
  done
else
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping fine mapping with CAVIAR"
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
fi
