#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Input arguments
chr    <- as.character(args[1])
phased <- as.character(args[2])
K      <- as.integer(args[3])

# Load libraries
source("../utils/util.R")
options(scipen = 999)

# Location of data
fp.dir <- "/scratch/PI/candes/ukbiobank_tmp/fastphase"

# Load true haplotypes from INP file
basename <- paste(fp.dir, "/test/ukb_hap_chr", chr, sep="")
Haplotypes.true <- read.inp(basename, progress=TRUE)
H.true <- as.matrix(select(Haplotypes.true, -Subject))

# Load masked haplotypes from INP file
basename <- paste(fp.dir, "/test/ukb_hap_chr", chr, "_masked", sep="")
Haplotypes.masked <- read.inp(basename, progress=TRUE)
H.masked <- as.matrix(Haplotypes.masked)

# Load imputed haplotypes from INP file
basename <- paste(fp.dir, "/test_imputed/", phased, "_K", K, "/ukb_hap_chr", chr, "_imputed", sep="")
Haplotypes.imputed <- read.inp(basename, progress=TRUE)
H.imputed <- as.matrix(Haplotypes.imputed)

# Compute imputation error
error <- mean(H.imputed[is.na(H.masked)]!=H.true[is.na(H.masked)])
cat(sprintf("Imputation error: %.4f\n", error))

# Save imputation error
Results = tibble("chr"=as.character(chr), "n"=as.integer(nrow(H.true)), "p"=as.integer(ncol(H.true)),
                 "miss"=round(mean(is.na(H.masked)),4), 
                 "phase"=phased, "K"=as.integer(K), "error"=round(error,4))

filename <- paste(fp.dir, "/test_imputed/", phased, "_K", K, "/ukb_hap_chr", chr, "_error.txt", sep="")
write_delim(Results, filename, delim="\t")
