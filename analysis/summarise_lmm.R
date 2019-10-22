#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

phenotype <- "bmi"
clumping <- "0.00000005"
fold <- "00"

args <- commandArgs(trailingOnly=TRUE)
phenotype <- as.character(args[1])
clumping <- as.character(args[2])
fold <- as.character(args[3])

library(tidyverse)
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
genotypes.dir <- "/scratch/PI/candes/ukbiobank/genotypes"
chr.list <- seq(1, 22)
resolution.dummy <- "Radj100"

# Load list of variants
Variants <- lapply(chr.list, function(chr) {
    cat(sprintf("Loading list of variants for chromosome %d...\n", chr))
    bim.file <- sprintf("%s/ukb_gen_chr%s.bim", genotypes.dir, chr)
    read_tsv(bim.file, col_names=c("CHR", "SNP", "X1", "BP", "A0", "A1"), col_types=cols())
})
Variants <- do.call("rbind", Variants)

# Load clumped results
if(fold=="00") {
    clumped.file <- sprintf("%s/analysis/bolt_clumped/%s_%s.tab", scratch, phenotype, clumping)
} else {
    clumped.file <- sprintf("%s/analysis/bolt_clumped/%s_fold_%s_%s.tab", scratch, phenotype, fold, clumping)
}
Clumped.raw <- read.table(clumped.file, header=TRUE) %>% as_tibble() %>%
    select(CHR, SNP, BP, SP2)

# Parse clumped results
if(nrow(Clumped.raw)>0) {
    # Extract lead SNPs
    Clumped.raw.lead <- Clumped.raw %>% mutate(SP2=SNP) %>% select(CHR, SNP, BP, SP2)
    # Extract secondary SNPs
    Clumped.raw.secondary <- Clumped.raw %>%
        mutate(SP2 = gsub("\\(1\\)", "", SP2)) %>%
        separate_rows(SP2, sep=",") %>%
        select(CHR, SNP, BP, SP2) %>%
        filter(SP2!="NONE")
    # Cross-reference with complete list of variants
    Clumped <- rbind(Clumped.raw.lead, Clumped.raw.secondary) %>%
        as_tibble() %>% arrange(CHR, BP, SP2) %>%
        mutate(SNP.lead=as.character(SNP), BP.lead=BP, SNP=as.character(SP2)) %>%
        select(CHR, SNP.lead, BP.lead, SNP) %>%
        left_join(select(Variants, CHR, SNP, BP), by=c('CHR', 'SNP')) %>%
        select(CHR, SNP.lead, BP.lead, SNP, BP)
} else {
    Clumped <- tibble()
}

# Load LMM p-values and add them to the list of clumped results
if(fold=="00") {
    lmm.file <- sprintf("%s/analysis/bolt/%s_stats.txt", scratch, phenotype)
} else {
    lmm.file <- sprintf("%s/analysis/bolt/%s_stats_fold_%s.txt", scratch, phenotype, fold)
}
LMM <- read_tsv(lmm.file, col_types=cols())
if(! "P" %in% colnames(LMM)) {
    if("P_BOLT_LMM" %in% colnames(LMM)) {
        LMM <- LMM %>% mutate(P=P_BOLT_LMM)
    } else {
        LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
    }
}
Clumped <- Clumped %>% inner_join(LMM %>% select(CHR,SNP,BP,P), by = c("CHR", "SNP", "BP"))

# Summarise results by clump
Discoveries <- Clumped %>% group_by(CHR, SNP.lead, BP.lead) %>%
    summarise(P=min(P), BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n())

# Save results
if(fold=="00") {
    regions.file <- sprintf("%s/discoveries/%s_lmm_regions.txt", scratch, phenotype)
} else {
    regions.file <- sprintf("%s/discoveries/%s_lmm_regions_fold_%s.txt", scratch, phenotype, fold)
}
Discoveries %>% write_delim(regions.file, delim=" ")
cat(sprintf("List of %d discovered regions written on:\n%s\n", nrow(Discoveries), regions.file))

if(fold=="00") {
    variants.file <- sprintf("%s/discoveries/%s_lmm_variants.txt", scratch, phenotype)
} else {
    variants.file <- sprintf("%s/discoveries/%s_lmm_variants_fold_%s.txt", scratch, phenotype, fold)
}
Clumped %>% write_delim(variants.file, delim=" ")
cat(sprintf("List of %d discovered variants written on:\n%s\n", nrow(Clumped), variants.file))
