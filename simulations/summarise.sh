#!/bin/bash

SUMMARY=$1
EXPERIMENT=$2
CONSERVATIVE=$3

mkdir -p "/scratch/PI/candes/ukbiobank_tmp/simulations/summary/"

source activate ukb
Rscript --vanilla "summarise_"$SUMMARY".R" $EXPERIMENT $CONSERVATIVE
