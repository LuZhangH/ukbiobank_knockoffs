#!/bin/bash

PHENOTYPE=$1
SEED=$2
CONSERVATIVE=$3

Rscript --vanilla multiresolution.R $PHENOTYPE $SEED $CONSERVATIVE
