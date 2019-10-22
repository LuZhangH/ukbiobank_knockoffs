#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Input arguments
chr <- args[1]

seed <- 2018

# Load libraries
suppressMessages(library(snpStats))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(SNPknock))
source("../utils/util.R")

# Location of holdout haplotypes
fp.dir <- "/scratch/PI/candes/ukbiobank_tmp/fastphase/test"

# Load haplotypes from INP file
basename <- paste(fp.dir, "/ukb_hap_chr", chr, sep="")
H <- read.inp(basename, progress=FALSE)
H <- as.matrix(select(H, -Subject))

# Generate mask and apply it
missingness <- 1/2
set.seed(seed)
n <- nrow(H)
p <- ncol(H)
mask <- matrix(rbinom(n*p,1,missingness),n,p)
mask.nnz <- which(mask!=0, arr.ind = T)
H[mask.nnz] <- NA

# Save the masked haplotypes
out.file <- paste(fp.dir, "/ukb_hap_chr", chr, "_masked.inp", sep="")
SNPknock.fp.writeX(H, phased=TRUE, out_file=out.file)
