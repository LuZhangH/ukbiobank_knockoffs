#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

args <- commandArgs(trailingOnly=TRUE)

# Default arguments
experiment <- "four"
n.signals <- 100
phenotype <- "Y_a10"
fold <- "01"
clumping <- "0.00000005"

# Input arguments
experiment <- as.character(args[1])
n.signals <- as.integer(args[2])
phenotype <- as.character(args[3])
fold <- as.character(args[4])
clumping <- as.character(args[5])

# Load libraries
suppressMessages(library(MASS))
source("utils/fine_mapping.R")
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))

K <- 50
resolution <- "Radj1"
chr.list <- seq(1, 22)

# Root directory
tmp.dir <- "/scratch/PI/candes/ukbiobank_tmp"

###########################
## Load list of variants ##
###########################

# Load list of variants
Variants <- lapply(chr.list, function(chr) {
    key.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.key", tmp.dir, resolution, K, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants.chr$Resolution <- resolution
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants)

# Load list of causal variants
causal.file <- sprintf("%s/simulations/phenotypes/%s_causal_s%d_%s.txt", tmp.dir, experiment, n.signals, fold)
Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, SNP, BP, Knockoff)
Variants <- Variants %>% mutate(Causal = (BP %in% Causal$BP))

Variants.original <- Variants %>% filter(Resolution==resolution, Knockoff==FALSE) %>%
    mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP))

Posterior <- lapply(chr.list, function(chr) {
    ##############################
    ## Load association results ##
    ##############################

    #################################
    ## Load list of clumps results ##
    #################################
    # Read clumping file
    clumped.file <- sprintf("%s/simulations/bolt/fold_%s/clumped/stats_%s_s100_%s_%s.clumped",
                          tmp.dir, fold, experiment, phenotype, clumping)
    Clumped.raw <- read.table(clumped.file, header=TRUE) %>% as_tibble() %>%
        select(CHR, SNP, BP, P, SP2) %>%
        filter(CHR==chr)

    if(nrow(Clumped.raw)==0) {
        return(tibble())
    }
    
    # Extract lead SNPs
    Clumped.raw.lead <- Clumped.raw %>% mutate(SP2=SNP) %>% select(CHR, SNP, BP, P, SP2)

    # Extract secondary SNPs
    Clumped.raw.secondary <- Clumped.raw %>%
        mutate(SP2 = gsub("\\(1\\)", "", SP2)) %>%
        separate_rows(SP2, sep=",") %>%
        select(CHR, SNP, BP, P, SP2) %>%
        filter(SP2!="NONE")

    # Cross-reference with complete list of variants
    Clumped <- rbind(Clumped.raw.lead, Clumped.raw.secondary) %>%
        as_tibble() %>% arrange(CHR, BP, SP2) %>%
        mutate(SNP.lead=as.character(SNP), BP.lead=BP, SNP=as.character(SP2)) %>%
        select(CHR, SNP.lead, BP.lead, SNP, P) %>%
        mutate(SNP.lead = gsub(".A", "", SNP.lead), SNP.lead = gsub(".B", "", SNP.lead)) %>%
        mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP)) %>%
        left_join(select(Variants.original, CHR, SNP, BP, Causal),
                  by=c('CHR', 'SNP')
                  ) %>%
        select(CHR, SNP.lead, BP.lead, SNP, BP, P, Causal)

        cat(sprintf("Loading marginal pvalues for chromosome %d ...\n", chr))
    
    # Read association file and remove knockoffs
    assoc.file <- sprintf("%s/simulations/pvalues/%s/fold_%s/s100_%s_chr%d.assoc.linear",
                          tmp.dir, experiment, fold, phenotype, chr)
    Association <- read.table(assoc.file, header=T) %>% as_tibble() %>%
        inner_join(select(Variants.original, CHR, SNP, BP, Causal), by = c("CHR", "SNP", "BP"))

    ###################
    ## Load LD table ##
    ###################
    cat(sprintf("Loading LD table for chromosome %d ...\n", chr))
    corr.filename <- sprintf("%s/stats/ukb_gen_chr%s.ld", tmp.dir, chr)
    LD.chr <- read_table2(corr.filename, col_types=cols())

    # Apply CAVIAR to each clump on this chromosome
    Posterior.chr <- fine.mapping.caviar(Clumped, Association, LD.chr) %>%
        left_join(select(Clumped, CHR, SNP, BP, Causal), by = c("CHR", "SNP", "BP"))

    return(Posterior.chr)
})
Posterior <- do.call("rbind", Posterior)

# Write posterior to file
out.file <- sprintf("%s/simulations/caviar/%s_s%d_%s_%s_%s.txt",
                    tmp.dir, experiment, n.signals, phenotype, fold, clumping)
Posterior %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))
