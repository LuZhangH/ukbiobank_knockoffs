#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Input arguments
chr <- args[1]
nsubjects <- args[2]

seed <- 123

# Load libraries
suppressMessages(library(snpStats))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
source("../utils/util.R")
source("util.R")

# Location of data
fp.dir <- paste("/scratch/PI/candes/ukbiobank_tmp/fastphase/CV_",nsubjects,sep="")

# Load haplotypes from INP file
basename <- paste(fp.dir, "/input/ukb_hap_chr", chr, "_test", sep="")
H <- read.inp(basename, progress=TRUE)
H <- as.matrix(select(H, -Subject))

# Generate mask
missingness <- 1/2
set.seed(seed)
n <- nrow(H)
p <- ncol(H)
mask <- matrix(rbinom(n*p,1,missingness),n,p)
mask.nnz <- which(mask!=0, arr.ind = T)

# Apply mask to haplotypes and save them in INP format
H[mask.nnz] <- NA
out.dir <- paste(fp.dir, "/input/mask01", sep="")
dir.create(out.dir, showWarnings = FALSE)
out.file <- paste(out.dir, "/ukb_hap_chr", chr, "_test.inp", sep="")
SNPknock.fp.writeX(H, phased=TRUE, out_file=out.file)
