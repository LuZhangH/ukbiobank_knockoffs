#!/bin/bash

EXPERIMENT_NAME=$1
EXPERIMENT_NUMB=$2
FOLD=$3
CONSERVATIVE=$4

Rscript --vanilla multiresolution.R $EXPERIMENT_NAME $EXPERIMENT_NUMB $FOLD $CONSERVATIVE
