#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Read input arguments
chr <- as.character(args[1])
clumping.factor <- as.numeric(args[2])
clumping.method <- as.character(args[3])
K <- as.numeric(args[4])
n.load <- as.integer(args[5])

## chr <- 22
## clumping.factor <- 10
## clumping.method <- "Radj"
## K <- 50
## n.load <- 1000

if(n.load == -1) n.load <- Inf

cat(paste("\nThis script is verifying that the knockoff-augmented genotypes for ",chr," have been packaged correctly.\n",
          "Running with K = ",K, ", clumping method ", clumping.method,
          " and clumping factor = ", clumping.factor/100, ".\n", sep=""))

cat("Loading R libraries... ")
suppressMessages(library(tidyverse))
suppressMessages(library(snpStats))
#suppressMessages(library(data.table))
source("../utils/util.R") # Load util functions
cat("done.\n")

# Specify location of data
pi.dir <- "/scratch/PI/candes"
rep.dir <- paste(pi.dir, "/ukbiobank_tmp/clumping/", clumping.method, clumping.factor, sep="")
fp.dir.in  <- paste(pi.dir, "/ukbiobank_tmp/fastphase/data", sep="")
knock.dir <- paste(pi.dir, "/ukbiobank_tmp/knockoffs/", clumping.method, clumping.factor, "_K", K, "/", sep="")

# Load augmented genotype matrix
basename <- paste(knock.dir, "ukb_gen_chr", chr, sep="")
Data <- read.plink(basename)
XA <- as(Data$genotypes, "numeric")

# Load knockoff swap key
key.file <- paste(basename, ".key", sep="")
Variants <- read_delim(key.file, delim=" ", col_types=cols(), progress=FALSE) %>%
    separate(Variant, c("Variant", "Suffix"), sep="\\.")

# Separate original and knockoff variables
idx.original <- which(Variants$Knockoff==FALSE)
idx.knockoff <- which(Variants$Knockoff==TRUE)
X  <- XA[,idx.original]
colnames(X) <- Variants$Variant[idx.original]
Xk <- XA[,idx.knockoff]
colnames(Xk) <- Variants$Variant[idx.knockoff]

# Define function to match reference alleles
match.alleles <- function(X, X.ref) {
    flipped <- sapply(1:ncol(X), function(j) {
        not.missing <- which(!is.na(X.ref[,j]))
        error.original <- sum(X[not.missing,j] != X.ref[not.missing,j])
        error.flipped  <- sum(X[not.missing,j] != 2-X.ref[not.missing,j])
        if(error.original == 0) return(0)
        else if(error.flipped == 0) return(1)
        else return(-1)
    })
    stopifnot(sum(flipped==-1)==0)
    cat(sprintf("%d alleles out of %d will be flipped. ", sum(flipped==1), length(flipped)))
    cat("No errors were found.\n")
    return(which(flipped==1))
}

# Load original genotypes
filename <- paste(fp.dir.in, "/ukb_gen_chr", chr, ".raw", sep="")
Genotypes.colnames <- strsplit(readLines(file(filename, open="r"), 1), "\t")[[1]]
Genotypes.colnames <- gsub("\\_.*", "", Genotypes.colnames)
Genotypes <- read_delim(filename, delim="\t", skip=1, n_max=n.load, progress=FALSE,
                        col_names=Genotypes.colnames,
                        col_types=cols(FID=col_character(), IID=col_character(),
                                       PAT=col_character(), MAT=col_character(),
                                       SEX=col_character(), PHENOTYPE=col_character(),
                                       .default = col_integer()))
X.o <- Genotypes %>% mutate(Subject=FID) %>% select(Variants$Variant[idx.original]) %>% as.matrix

# Verify compatibility of genotypes
cols.flip <- match.alleles(X, X.o)
