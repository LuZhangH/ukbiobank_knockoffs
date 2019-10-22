#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

# Install packages bigsnpr and bigstatsr
# devtools::install_github("privefl/bigstatsr")
# devtools::install_github("privefl/bigsnpr")

# Documentation here:
# https://privefl.github.io/bigsnpr/reference/index.html
# https://privefl.github.io/bigstatsr/reference/big_spLinReg.html

# Load packages
library(tidyverse)
library(devtools)
library(bigsnpr)

# Default arguments
resolution <- "Radj100"
seed <- 1

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
resolution <- as.character(args[1])
seed <- as.character(args[2])

# Other parameters
scratch <- "/scratch/PI/candes/ukbiobank_tmp"

####################
## Load genotypes ##
####################

tmp.file <- sprintf("%s/augmented_data_big/ukb_gen_%s_s%s", scratch, resolution, seed)
rds.file <- sprintf("%s.rds", tmp.file)

if(file.exists(rds.file)){
    cat(sprintf("Found FBM in %s. Skipping. \n", rds.file))
} else {
    bed.file <- sprintf("%s/augmented_data_big/ukb_gen_%s_s%s.bed", scratch, resolution, seed)
    cat(sprintf("Making temporary file for FBM: %s\n", sprintf("%s.bk", tmp.file)))
    cat("Reading genotype file and creating FBM object... ")
    rds <- snp_readBed(bed.file, backingfile = tmp.file)
    cat("done.\n")
}

# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
obj.bigSNP <- snp_attach(rds.file)
cat("done.\n")
