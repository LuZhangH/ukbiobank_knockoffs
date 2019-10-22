#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Read input arguments
chr <- as.character(args[1])
clumping.factor <- as.integer(args[2])
clumping.method <- as.character(args[2])
K <- as.numeric(args[4])

# Global parameters
n.load <- -1
n.skip <- 0
block.size <- 1000      # Number of individuals to process simultaneously
phased.flag <- "phased"

cat(paste("This script is generating Gaussian knockoffs for chromosome ",chr,".\n",
          "Running with K = ",K, ", clumping method ", clumping.method,
          " and clumping factor = ", clumping.factor/100, ".\n", sep=""))

cat("Loading R libraries... ")
suppressMessages(library(tidyverse))
suppressMessages(library(SNPknock))
suppressMessages(library(doParallel))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(knockoff))
source("../utils/util.R") # Load util functions
cat("done.\n")

# Specify location of data
pi.dir <- "/scratch/PI/candes"
clump.dir <- paste(pi.dir, "/ukbiobank_tmp/clumping", sep="")
rep.dir <- paste(clump.dir, "/", clumping.method, clumping.factor, sep="")
knock.dir <- paste(pi.dir, "/ukbiobank_tmp/knockoffs_gau/", clumping.method, clumping.factor, "/", sep="")
dir.create(knock.dir, showWarnings = FALSE)

# Load list of prototype variants
cat("Loading list of variant prototypes... ")
rep.file <- paste(rep.dir, "/rep_chr", chr, ".txt", sep="")
Variants <- read_delim(rep.file, delim=" ", progress = FALSE, col_types = cols())
cat("done.\n")

# Load correlation matrix computed with PLINK
cat("Loading correlation matrix computed with PLINK... ")
filename <- paste(pi.dir, "/ukbiobank_tmp/stats/ukb_gen_chr", chr, ".ld", sep="")
LD <- read_table(filename, col_types=cols(CHR_A = col_integer(), BP_A = col_integer(),
                                          SNP_A = col_character(), CHR_B = col_integer(),
                                          BP_B = col_integer(), SNP_B = col_character(),
                                          R2 = col_double(), DP = col_double()))
LD <- LD %>% filter(BP_A %in% Variants$Position, BP_B %in% Variants$Position) %>% mutate(R=sqrt(R2))
Sigma <- ld.to.mat(LD, clumping.method, Variants$Position)
cat("done.\n")

# Prepare to generate Gaussian knockoffs
cat("Solving SDP for Gaussian knockoffs... \n")
diag_s = create.solve_asdp(as.matrix(Sigma), max.size=500, verbose=TRUE)
cat("done.\n")
# Save diag_s to file
knock.file <- paste(knock.dir, "ukb_gen_chr", chr, ".diags", sep="")
write(diag_s, ncolumns=1, knock.file)
