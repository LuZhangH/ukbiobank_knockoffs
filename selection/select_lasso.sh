#!/bin/bash
# UK Biobank GWAS
#
# Class: script
# 
# Description 
#
# Authors: Matteo Sesia
# Date:  Sept 24, 2018

# Input parameters
PHENOTYPE=$1
CLUMPING=$2

# Load libraries
ml gcc
ml R/3.5.1

# Call R script
Rscript --vanilla select_lasso.R $PHENOTYPE $CLUMPING
