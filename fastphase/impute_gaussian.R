#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Input arguments
chr    <- as.character(args[1])
phased <- as.character(args[2])
K      <- as.integer(args[3])
ld.measure <- "R"

# Load libraries
library(Matrix)
library(corpcor)
source("../utils/util.R")
options(scipen = 999)

# Location of data
fp.dir <- "/scratch/PI/candes/ukbiobank_tmp/fastphase"
stats.dir <-"/scratch/PI/candes/ukbiobank_tmp/stats"
    
# Load true haplotypes from INP file
basename <- paste(fp.dir, "/test/ukb_hap_chr", chr, sep="")
Haplotypes.true <- read.inp(basename, progress=TRUE)
H.true <- as.matrix(select(Haplotypes.true, -Subject))

# Load list of variants
Variants <- read_delim(paste(basename, ".legend", sep=""), delim = " ", skip=1,
                           col_names=c("Variant", "Position", "A1", "A2"))
frq.filename = sprintf("%s/stats/ukb_gen_chr%s.frq", tmp.dir, chr)

# Load sparse covariance matrix
corr.filename = sprintf("%s/ukb_gen_chr%s.ld", stats.dir, chr)
LD = read_table(corr.filename)
LD = filter(LD, SNP_A %in% Variants$Variant, SNP_B %in% Variants$Variant)

# Load masked haplotypes from INP file
basename <- paste(fp.dir, "/test/ukb_hap_chr", chr, "_masked", sep="")
Haplotypes.masked <- read.inp(basename, progress=TRUE)
H.masked <- as.matrix(Haplotypes.masked)

impute.block <- function(H.true.block, H.masked.block, Variants.block) {
    # Convert LD table to a symmetric square matrix
    LD.block = filter(LD, SNP_A %in% Variants.block$Variant, SNP_B %in% Variants.block$Variant)
    Sigma <- ld.to.mat(LD.block, ld.measure, Variants.block$Position)

    # Convert correlation matrix to covariance matrix
    stdevs <- sqrt(diag(cov(H.true.block)))
    Sigma <- stdevs %*% t(stdevs) * Sigma

    # Compute marginal expectations
    mu <- colMeans(H.true.block)

    # Impute missing values
    H.imputed.block <- matrix(NA, nrow(H.masked.block), ncol(H.masked.block))
    for(i in 1:nrow(H.true.block)) {
        missing <- which(is.na(H.masked.block[i,]))
        mu.cond <- mu[missing] + (H.masked.block[i,-missing] - mu[-missing]) %*%
            pseudoinverse(Sigma[-missing,-missing]) %*% Sigma[-missing,missing]
        H.imputed.block[i,missing] <- as.numeric(round(mu.cond))
        H.imputed.block[i,-missing] <- H.masked.block[i,-missing]
    }

    # Compute imputation error
    error <- mean(H.imputed.block[is.na(H.masked.block)]!=H.true.block[is.na(H.masked.block)])
    return(error)
}

# Divide variants into smaller blocks
block.size <- 100
chunk2 <- function(x,n) split(x, cut(seq_along(x), round(length(x)/n), labels = FALSE)) 
blocks <- chunk2(1:nrow(Variants), block.size)

# Impute missing values, block-by-block
imputation.errors <- rep(NA, length(blocks))
pb <- txtProgressBar(1, length(blocks), style=3)
setTxtProgressBar(pb, 0)
for(b in 1:length(blocks)) {
    block <- blocks[[b]]
    H.true.block <- H.true[,block]
    H.masked.block <- H.masked[,block]
    Variants.block <- Variants[block,]
    imputation.errors[b] <- impute.block(H.true.block, H.masked.block, Variants.block)
    setTxtProgressBar(pb, b)
}

# Save results
Results <- tibble(n=c(nrow(H.true)), K="Gaussian", error=round(median(imputation.errors),5))
out.file <- paste(fp.dir, "/test_imputed/", phased, "_gaussian/ukb_hap_chr", chr, "_error.txt", sep="")
write_delim(Results, out.file, delim="\t")

