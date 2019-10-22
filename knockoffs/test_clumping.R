#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

chr = 22
clumping.method <- "Radj"
    
suppressMessages(library(adjclust))
library(tidyverse)
suppressMessages(library(Matrix))
source("../utils/util.R")

# Parse input
ld.measure <- substring(clumping.method, 1, 1)
clustering <- substring(clumping.method, 2)

# Paths to data
dat.dir <- "/scratch/PI/candes/ukbiobank"
tmp.dir <- "/scratch/PI/candes/ukbiobank_tmp"

# Load list of variants for phased haplotypes that have passed QC
positions.filename <- paste(tmp.dir, "/fastphase/data/ukb_hap_chr", chr, ".legend", sep="")
variant.names <- read_delim(positions.filename, delim=" ", skip=1, col_types=cols(),
                           col_names=c("Variant", "Position", "A1", "A2"))
bim.filename <- paste(dat.dir, "/genotypes/ukb_gen_chr", chr, ".bim", sep="")
Variants <- as_tibble(read.table(bim.filename, sep="", header=FALSE, stringsAsFactors=FALSE,
                                col.names=c("Chr", "Variant", "X1", "Position", "A1", "A2")))
Variants <- semi_join(Variants, variant.names, by=c("Variant"))

# Load MAF
frq.filename <- sprintf("%s/stats/ukb_gen_chr%s.frq", tmp.dir, chr)
Frequencies <- as_tibble(read.table(frq.filename, sep="", header=TRUE, stringsAsFactors=FALSE,
                                   col.names=c("Chr","Variant","A1","A2","MAF","NCHROBS")))
Variants <- inner_join(Variants, select(Frequencies, c("Variant", "MAF", "NCHROBS")), by=c("Variant"))

# Compute dendrogram if not already present
dend.filename <- sprintf("%s/clumping/%s/grp_chr%s.RData", tmp.dir, clumping.method, chr)
if(file.exists(dend.filename)) {
    cat(sprintf("File %s found.\n", dend.filename))
    cat("Loading clustering dendrogram... ")
    load(dend.filename)
    cat("done.\n")
}

# Choose groups by cutting the dendrogram
resolution.list <- c(1,2,5,10,20,50)/100
groups.list <- sapply(resolution.list, function(resolution) {
    cutree(Sigma.clust, k = round(resolution*nrow(Variants)))
})
colnames(groups.list) <- resolution.list

#groups.list <- groups.list[1:1000,]
all.contained <- sapply(seq(ncol(groups.list)-1), function(j) {
    is.contained <- sapply(unique(groups.list[,j+1]), function(val) {
        val.idx <- which(groups.list[,j+1]==val)
        contained <- length(unique(groups.list[val.idx,j]!=val))==1
    })
    sum(!is.contained)==0
})

sum(!all.contained)==0
