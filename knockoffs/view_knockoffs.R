#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(snpStats))
suppressMessages(library(tidyverse))
source("../utils/util.R")

# Specify chromosome
chr <- 22
cutoff <- 50
K <- 30

n.load <- 1000
n.skip <- 0

# Specify location of data
pi.dir <- "/scratch/PI/candes"
rep.dir <- paste(pi.dir, "/ukbiobank_tmp/clumping/R", cutoff, sep="")
data.dir <- paste(pi.dir,"/ukbiobank/genotypes/", sep="")
fp.dir <- paste(pi.dir, "/ukbiobank_tmp/fastphase/phased_K", K, "/", sep="")
knock.dir <- paste(pi.dir, "/ukbiobank_tmp/knockoffs/R", cutoff, "_K", K, "/", sep="")

# Load list of prototype variants
cat("Loading list of variant prototypes... ")
rep.file <- paste(rep.dir, "/rep_chr", chr, ".txt", sep="")
Variants <- read_delim(rep.file, delim=" ", progress = FALSE, col_types = cols())
ref.file <- paste(pi.dir, "/ukbiobank_tmp/fastphase/data/ukb_hap_chr", chr, ".ref", sep="")
Variants.ref <- read_delim(ref.file, delim="\t", progress = FALSE, col_types = cols())
Variants <- left_join(Variants, Variants.ref, by="Variant")
cat("done.\n")

# Load list of subjects
cat("Loading list of subjects... ")
sam.file <- paste(knock.dir, "ukb_hap_chr",chr,".sample", sep = "")
Subjects <- read_delim(sam.file, delim=" ", skip=n.skip+1, n_max=n.load, progress = TRUE,
                       col_types=cols(.default = col_integer()), col_names=c("pedigree"))
cat("done.\n")

# Load genotype data
cat("Loading genotypes... ")
basename <- paste(pi.dir, "/ukbiobank_tmp/fastphase/data/ukb_gen_chr", chr, sep="")
Genotypes <- read.raw(basename, n.skip=n.skip, n.load=n.load, progress=TRUE)
subjects.gen <- Genotypes$Subject
Genotypes <- select(Genotypes, Variants$Variant)
missing.gen <- mean(is.na(Genotypes))
X <- as.matrix(Genotypes)
cat("done.\n")
# Impute missing genotypes
for(j in 1:ncol(X)) {
    X[is.na(X[,j]),j] = median(X[,j], na.rm=TRUE)
}

# Load knockoffs
cat("Loading knockoff genotypes... ")
knock.file <- paste(knock.dir, "ukb_hap_chr", chr, ".haps.t", sep="")
Hk.dat <- read_delim(knock.file, delim=" ", skip=n.skip, n_max=2*n.load, col_names=FALSE,
                   col_types=cols(.default = col_integer()),progress=TRUE)
Hk <- as.matrix(Hk.dat)
Xk <- Hk[seq(1,nrow(Hk),by=2),] + Hk[seq(2,nrow(Hk),by=2),]
cat("done.\n")

# Match reference alleles
flip.sign.X <- which(Variants$Flip==1)
Xk[,flip.sign.X] <- 2 - Xk[,flip.sign.X]

# Compare SNP column means
plot(colMeans(X),colMeans(Xk),col=rgb(0,0,0,alpha=0.1), pch=16,cex=1); abline(a=0, b=1, col='red', lty=2)

# Compare correlations between consecutive SNPs
corrX = sapply(2:dim(X)[2], function(j) cor(X[,j-1],X[,j]))
corrXk = sapply(2:dim(X)[2], function(j) cor(Xk[,j-1],Xk[,j]))
plot(abs(corrX),abs(corrXk),col=rgb(0,0,0,alpha=0.1), pch=16,cex=1); abline(a=0, b=1, col='red', lty=2)

# Compare correlations between original SNPs and their successive knockoff
corrXXk = sapply(2:dim(X)[2], function(j) cor(X[,j-1],Xk[,j]))
plot(corrX,corrXXk,col=rgb(0,0,0,alpha=0.1), pch=16,cex=1); abline(a=0, b=1, col='red', lty=2)

# Measure knockoff quality
knock.quality = sapply(1:dim(X)[2], function(j) cor(X[,j],Xk[,j]))
hist(knock.quality)
