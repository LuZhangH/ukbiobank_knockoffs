#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

args <- commandArgs(trailingOnly=TRUE)

# Load libraries
suppressMessages(library(MASS))
source("../simulations/utils/fine_mapping.R")
source("../utils/utils_clumping.R")
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))

# Default arguments
experiment.name <- "six"
amplitude <- 10
fold <- 1
clumping <- "0.00000005"

# Input arguments
experiment.name <- as.character(args[1])
amplitude <- as.character(args[2])
fold <- as.character(args[3])
clumping <- as.character(args[4])

# Other parameters
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1,22)
resolution <- "Radj100"
resolution.list <- c(resolution)
pheno.name <- sprintf("Y_a%s_%s", amplitude, fold)

###########################
## Load list of variants ##
###########################

# Load list of variants
Variants <- lapply(chr.list, function(chr) {
    key.file <- sprintf("%s/knockoffs/%s_K50/ukb_gen_chr%d.key", scratch, resolution, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants.chr$Resolution <- resolution
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants)

# Load list of causal variants
causal.file <- sprintf("%s/simulations_small/phenotypes/%s_causal.txt", scratch, experiment.name)
Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, BP, Locus, Sign, Scale)
Variants <- Variants %>%
    left_join(mutate(Causal,Causal=TRUE),by = c("CHR", "BP")) %>%
    mutate(Causal=replace_na(Causal,FALSE)) %>%
    mutate(Causal=ifelse(Knockoff, FALSE, Causal))

# Extract list of original variants
Variants.original <- Variants %>%
    filter(Resolution==resolution, Knockoff==FALSE) %>%
    mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP))
                                       
##################
## Fine mapping ##
##################

Posterior <- lapply(chr.list, function(chr) {
    ##############################
    ## Load association results ##
    ##############################

    #################################
    ## Load list of clumps results ##
    #################################
    # Read clumping file
    clumped.file <- sprintf("%s/simulations_small/bolt_clumped/fold_%s/stats_%s_%s_%s.clumped",
                          scratch, fold, experiment.name, amplitude, clumping)
    if(file.exists(clumped.file)) {
        Clumped.raw <- read.table(clumped.file, header=TRUE) %>% as_tibble() %>%
            select(CHR, SNP, BP, P, SP2) %>%
            filter(CHR==chr)

        # Parse list of clumped discoveries
        if(nrow(Clumped.raw)==0) {
            Clumped <- tibble()
        } else {
            Variants.prefix <- Variants.original %>%
                mutate(SNP.sec = gsub(".A", "", SNP), SNP.sec = gsub(".B", "", SNP.sec), BP.sec=BP) %>%
                select(CHR, SNP, SNP.sec, BP.sec)
            # Extract lead SNPs
            Clumped.raw.lead <- Clumped.raw %>% select(CHR, BP, P) %>%
                inner_join(select(Variants.original, CHR, SNP, BP), by = c("CHR", "BP")) %>%
                mutate(SNP.lead=SNP, BP.lead=BP, SNP=SNP) %>%
                select(CHR, SNP.lead, BP.lead, P, SNP, BP)
            # Extract secondary SNPs
            Clumped.raw.secondary <- Clumped.raw %>%
                mutate(SNP.lead=SNP, SNP.sec = gsub("\\(1\\)", "", SP2)) %>%
                separate_rows(SNP.sec, sep=",") %>%
                select(CHR, BP, P, SNP.sec) %>%
                filter(SNP.sec!="NONE") %>%
                inner_join(select(Variants.original, CHR, SNP, BP), by = c("CHR", "BP")) %>%
                mutate(SNP.lead=SNP, BP.lead=BP) %>%
                select(CHR, SNP.lead, BP.lead, P, SNP.sec) %>%
                inner_join(Variants.prefix, by = c("CHR", "SNP.sec")) %>%
                mutate(BP=BP.sec) %>%
                select(CHR, SNP.lead, BP.lead, P, SNP, BP, P)
            # Combine lists of lead and secondary SNPs
            Clumped <- rbind(Clumped.raw.lead, Clumped.raw.secondary) %>% as_tibble() %>%
                arrange(CHR, BP.lead, BP)
            # Cross-reference with complete list of variants
            Clumped <- Clumped %>%
                inner_join(select(Variants.original,CHR,SNP,BP,Causal), by = c("CHR", "SNP", "BP"))
        }
    } else {
        Clumped <- tibble()
    }

    # Consolidate clumps
    if(nrow(Clumped)>0) {
        Clumped.consolidated <- Clumped %>%
            group_by(CHR, SNP.lead, BP.lead) %>%
            summarise(P=min(P), BP.min=min(BP), BP.max=max(BP), Size=n(), Causal=any(Causal)) %>%
            mutate(BP.width=BP.max-BP.min) %>%
            ungroup() %>%
            consolidate_clumps(gap=1e5)
        Clumped.consolidated <- full_join(select(Clumped, CHR, SNP, BP, P, Causal),
                                          select(Clumped.consolidated, CHR, SNP.lead, BP.lead, BP.min, BP.max),
                                          by = c("CHR")) %>%
            filter(BP>=BP.min, BP<=BP.max)
    } else {
        Clumped.consolidated <- tibble()
    }

    if(nrow(Clumped.consolidated)>0) {
        # Read association file and remove knockoffs
        assoc.file <- sprintf("%s/simulations_small/pvalues/fold_%s/%s_%s_chr%d.assoc.linear",
                              scratch, fold, experiment.name, amplitude, chr)
        Association <- read.table(assoc.file, header=T) %>% as_tibble() %>%
            mutate(SNP=as.character(SNP)) %>%
            inner_join(select(Variants.original, CHR, SNP, BP, Causal), by = c("CHR", "SNP", "BP"))

        ###################
        ## Load LD table ##
        ###################
        cat(sprintf("Loading LD table for chromosome %d ...\n", chr))
        corr.filename <- sprintf("%s/stats/ukb_gen_chr%s.ld", scratch, chr)
        LD.chr <- read_table2(corr.filename, col_types=cols())

        # Apply CAVIAR to each clump on this chromosome
        Posterior.chr <- fine.mapping.caviar(Clumped.consolidated, Association, LD.chr) %>%
            left_join(select(Clumped.consolidated, CHR, SNP, BP, Causal), by = c("CHR", "SNP", "BP"))

    } else {
        Posterior.chr <- tibble()
    }
    
    return(Posterior.chr)
})
Posterior <- do.call("rbind", Posterior)

# Write posterior to file
out.file <- sprintf("%s/simulations_small/caviar/fold_%s/%s_%s_%s.txt",
                    scratch, fold, experiment.name, amplitude, clumping)
Posterior %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))
