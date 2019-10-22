#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")
suppressMessages(library(data.table))
suppressMessages(library(bigstatsr))

## install_bitbucket("msesia/hmm_knockoffs", ref = "master", subdir = "/SNPknock_R/SNPknockG",, auth_user = "msesia", password = "XXX")
## devtools::install("/home/users/msesia/Workspace/hmm_knockoffs/SNPknock_R/SNPknockG", build=T)

## Default parameters
chr <- 22

## Input parameters
args = commandArgs(trailingOnly=TRUE)
chr <- as.character(args[1])

## Data location
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
#scratch <- "/home/groups/candes/scratch/ukbiobank_tmp"

hapst.file <- sprintf("%s/fastphase/data/ukb_hap_chr%s.haps.t", scratch, chr)
bk.file <- sprintf("%s/knockoffs/data/ukb_hap_chr%s", scratch, chr)

## Load HAPS.T file
cat(sprintf("Loading transposed haps file...\n"))
H <- fread(hapst.file, sep=" ", showProgress=TRUE)
cat("Done.\n")

## Convert to FBM
cat(sprintf("Converting to FBM... "))
fbm <- as_FBM(H, type="unsigned short", backingfile=bk.file)
fbm$save()
cat("done.\n")
cat(sprintf("Written FBM: %s.rds\n", bk.file))
