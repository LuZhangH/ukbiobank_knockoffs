#!/bin/bash
#
# Dowload Po-Ru's p-values
#

# What should we do?
DO_DOWNLOAD=0

# Main storage location
META_DIR="/scratch/groups/candes/ukbiobank_tmp/meta"
mkdir -p $META_DIR
mkdir -p $META_DIR"/stats/"

PHENOTYPE_LIST=("body_HEIGHTz" "body_BMIz" "blood_PLATELET_COUNT" "bp_SYSTOLICadjMEDz" "disease_CARDIOVASCULAR" "disease_HYPOTHYROIDISM_SELF_REP")

# Remote server
REMOTE="https://data.broadinstitute.org/alkesgroup/UKBB/"

for PHENOTYPE in "${PHENOTYPE_LIST[@]}"; do    
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Downloading data for "$PHENOTYPE
  echo "----------------------------------------------------------------------------------------------------"
  echo ""
  
  # Dowload archive
  PACKAGE=$PHENOTYPE".sumstats"
  if [[ ! -f $META_DIR"/stats/"$PACKAGE ]]; then
    echo "Downloading package '$PACKAGE'.gz from remote sever ..."
    wget $REMOTE$PACKAGE".gz"
  else
    echo "Found "$META_DIR"/stats/"$PACKAGE
    echo "Skipping."
    continue
  fi

  # Unpack archive
  if [[ -f $PACKAGE".gz" ]]; then
    echo "Unpacking '$PACKAGE' ..."
    gunzip $PACKAGE".gz"
  else
    if [[ -f $PACKAGE ]]; then
      echo "Found "$PACKAGE
    else
      echo "Could not find "$PACKAGE
      echo "Skipping."
      continue
    fi        
  fi

  # Move archive to scratch partition
  if [[ -f $PACKAGE ]]; then
    echo "Moving '$PACKAGE' to: "$META_DIR"/stats/ ..."
    mv $PACKAGE $META_DIR"/stats/"
  fi

done
