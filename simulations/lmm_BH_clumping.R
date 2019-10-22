#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

# Main parameters (default)
experiment.name <- "seven"
amplitude <- 10

# Main parameters (input)
args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])
amplitude <- as.numeric(args[2])

# Load libraries
suppressMessages(library(tidyverse))

K <- 50
resolution <- "Radj1"
chr.list <- seq(1, 22)
fold <- "01"

# Root directory
scratch <- "/scratch/PI/candes/ukbiobank_tmp"

###########################
## Load list of variants ##
###########################

# Load list of variants
Variants <- lapply(chr.list, function(chr) {
    key.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.key", scratch, resolution, K, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants.chr$Resolution <- resolution
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants)

Variants.original <- Variants %>% filter(Resolution==resolution, Knockoff==FALSE) %>%
    mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP))

#######################
## Load LMM p-values ##
#######################

lmm.file <- sprintf("%s/simulations/bolt/fold_1/stats_%s_%s.txt", scratch, experiment.name, amplitude)
LMM <- read_tsv(lmm.file, col_types=cols())
if(! "P" %in% colnames(LMM)) {
    if("P_BOLT_LMM" %in% colnames(LMM)) {
        LMM <- LMM %>% mutate(P=P_BOLT_LMM)
    } else {
        LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
    }
}

#########################################
## Compute rejection threshold with BH ##
#########################################

BH.threshold <- function(pvalues, alpha) {
    m <- length(pvalues)
    pvalues.s <- sort(pvalues)
    ref.line <- (alpha/m)*(1:m)
    is.below <- which(pvalues.s<=ref.line)
    if(length(is.below)==0) {
        threshold <- 5e-8
    } else {
        threshold <- pvalues.s[max(is.below)]
    }
    return(threshold)
}

alpha <- 0.1
fdr.threshold <- BH.threshold(LMM$P, alpha)
cat(sprintf("The BH threshold at level %.2f is %f\n", alpha, fdr.threshold))

#############################
## Run clumping with PLINK ##
#############################

for(chr in chr.list) {

    plink.bfile <- sprintf("/scratch/PI/candes/ukbiobank/genotypes/ukb_gen_chr%d", chr)
    plink.exclude <- "/scratch/PI/candes/ukbiobank_tmp/simulations/variants/ukb_gen_exclude.txt"
    plink.stats <- sprintf("%s/simulations/bolt/fold_1/stats_%s_%s.txt.tmp",
                           scratch, experiment.name, amplitude)
    plink.threshold <- sprintf("%s", fdr.threshold)
    clump.dir <- sprintf("%s/simulations/bolt_clumped/fold_1/", scratch)
    plink.out <- sprintf("%s/stats_%s_%s_BH_chr%d", clump.dir, experiment.name, amplitude, chr)

    plink.command <- sprintf("plink --bfile %s \ --exclude %s \ --clump %s \ --clump-p1 %s \ --clump-r2 0.01 \ --clump-kb 5000 \ --out %s", plink.bfile, plink.exclude, plink.stats, plink.threshold, plink.out)

    if(file.exists(sprintf("%s.clumped", plink.out))) {
        cat(sprintf("Skipping clumping for chromosome %d.\n", chr))
    } else {
        system(plink.command)
        plink.clean.command <- sprintf("rm %s.log %s.nosex", plink.out, plink.out) 
        system(plink.clean.command)
    }
    
}
