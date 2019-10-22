#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(snpStats)
library(tidyverse)
source("util.R")
source("../utils/util.R")

# Specify chromosome and cutoff
chr <- 22
cutoff <- 50

n.load <- 1000
n.skip <- 0
K <- 20
ld.measure="R"

# Specify location of data
pi.dir <- "/scratch/PI/candes"
rep.dir <- paste(pi.dir, "/ukbiobank_tmp/clumping/", ld.measure, cutoff, sep="")
fp.dir <- paste(pi.dir, "/ukbiobank_tmp/fastphase/",sep="")

# Load haplotypes
basename <- paste(fp.dir, "data/ukb_hap_chr", chr, sep="")
Haplotypes <- read.inp(basename, n.skip=n.skip, n.load=n.load, progress=TRUE)
H.hap <- as.matrix(select(Haplotypes, -Subject))
X.hap <- H.hap[seq(1,nrow(H.hap),by=2),] + H.hap[seq(2,nrow(H.hap),by=2),]

# Load genotypes
basename <- paste(fp.dir, "data/ukb_gen_chr", chr, sep="")
Genotypes <- read.raw(basename, n.skip=n.skip, n.load=n.load, progress=TRUE)
X <- as.matrix(select(Genotypes, colnames(Haplotypes)[-1]))

flip.file <- paste(fp.dir, "data/ukb_hap_chr", chr, ".ref", sep="")
Flipped <- read_delim(flip.file, delim="\t")
flip.sign.X <- which(Flipped$Flip==1)
X.hap[,flip.sign.X] <- 2-X.hap[,flip.sign.X]

# Compare genotypes and haplotypes
plot(colMeans(X),colMeans(X.hap))
