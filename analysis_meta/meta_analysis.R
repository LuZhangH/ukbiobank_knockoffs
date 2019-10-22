#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")
library(tidyverse)
source("../utils/utils_clumping.R")

phenotype.code <- "body_HEIGHTz"
use.small <- TRUE
use.BH <- TRUE
    
args <- commandArgs(trailingOnly=TRUE)
phenotype.code <- as.character(args[1])
use.small <- as.logical(as.numeric(args[2]))
use.BH <- as.logical(as.numeric(args[3]))

if(use.BH && !use.small) {
    cat("BH results are not available for large analysis.\n")
    quit()
}

significance.meta <- "0.000000005"
resolution <- "Radj2"

# Define phenotype table
pheno.codes <- c("body_HEIGHTz",
                 "body_BMIz",
                 "blood_PLATELET_COUNT",
                 "bp_SYSTOLICadjMEDz",
                 "disease_CARDIOVASCULAR",
                 "disease_HYPOTHYROIDISM_SELF_REP"
                 )
pheno.names <- c("height", "bmi", "platelet", "sbp", "cvd", "hypothyroidism")
Phenotype.table <- tibble(Code=pheno.codes, Name=pheno.names)

# Match phenotype code
phenotype <- filter(Phenotype.table, Code==phenotype.code)$Name

scratch <- "/scratch/PI/candes/ukbiobank_tmp"
scratch.out <- sprintf("%s/meta/summary", scratch)
#scratch.out <- "/home/users/msesia/Workspace/ukbiobank/analysis_meta/results"

###################################
## Load results of meta analysis ##
###################################

meta.file <- sprintf("%s/meta/stats/%s.sumstats", scratch, phenotype.code)
Meta.pvalues <- read_tsv(meta.file, col_types=cols())
Meta.pvalues <- Meta.pvalues %>% mutate(BP=POS) %>% select(CHR, SNP, BP, P)

##################################
## Load results of KnockoffZoom ##
##################################
if(use.small) {
    discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_%s_fold_01.txt", scratch, phenotype, resolution)
} else {
    discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_%s.txt", scratch, phenotype, resolution)
}
Knockoffs <- read_delim(discoveries.file, delim=" ", col_types=cols())

###################################
## Knockoff discoveries vs. meta ##
###################################
match_to_meta <- function(chr, bp.min, bp.max, Meta.pvalues, significance) {
    gap <- 1e5
    matches <- Meta.pvalues %>%
        filter(CHR==chr, BP<=bp.max+gap, BP>=bp.min-gap, P<=significance) %>%
        arrange(P)
    if(nrow(matches)>0) {
        return(matches$BP[1])
    } else {
        return(NA)
    }
}
verify_with_meta <- function(Discoveries, Meta.pvalues, significance=5e-8) {
    # Remove rows that are not significant from results of meta analysis
    Meta.pvalues <- Meta.pvalues %>% filter(P<=significance)
    # Match each row of the discoveries with the meta analysis
    Discoveries <- Discoveries %>%
        rowwise() %>%
        mutate(Verified.BP = match_to_meta(CHR, BP.min, BP.max, Meta.pvalues, significance)) %>%
        ungroup()
    # Add SNP info
    Meta.pvalues.ref <- Meta.pvalues %>%
        mutate(Verified.BP=BP, Verified.SNP=SNP)
    Discoveries <- Discoveries %>%
        left_join(Meta.pvalues.ref, by=c("CHR", "Verified.BP"))
    return(Discoveries)
}

Knockoffs.confirmed <- Knockoffs %>% verify_with_meta(Meta.pvalues, significance=parse_number(significance.meta))

cat(sprintf("Verified knockoff discoveries: %d out of %d (%.2f%%).\n",
            sum(!is.na(Knockoffs.confirmed$Verified.BP)),
            nrow(Knockoffs.confirmed),
            100*mean(!is.na(Knockoffs.confirmed$Verified.BP))))

# Save results
if(use.small) {
    res.file <- sprintf("%s/knockoffs_%s_fold_01.txt", scratch.out, phenotype)
} else {
    res.file <- sprintf("%s/knockoffs_%s.txt", scratch.out, phenotype)
}
Knockoffs.confirmed %>% write_delim(res.file, delim=" ")
cat(sprintf("List of verified knockoff discoveries written on:\n %s\n", res.file))

##############################
## Load our LMM discoveries ##
##############################

if(use.small) {
    # Load clumped discoveries
    if(use.BH) {
        clumps.file <- sprintf("%s/discoveries/%s_lmm_regions_fold_01_BH.txt", scratch, phenotype)
    } else {
        clumps.file <- sprintf("%s/discoveries/%s_lmm_regions_fold_01.txt", scratch, phenotype)
    }
    LMM <- read_delim(clumps.file, delim=" ", col_types=cols())
    # Try to replicate
    LMM.confirmed <- LMM %>% verify_with_meta(Meta.pvalues, significance=parse_number(significance.meta))
    cat(sprintf("Verified LMM discoveries: %d out of %d (%.2f%%).\n",
            sum(!is.na(LMM.confirmed$Verified.BP)),
            nrow(LMM.confirmed),
            100*mean(!is.na(LMM.confirmed$Verified.BP))))
    # Save results
    if(use.BH) {
        res.file <- sprintf("%s/lmm_%s_fold_01_BH.txt", scratch.out, phenotype)
    } else {
        res.file <- sprintf("%s/lmm_%s_fold_01.txt", scratch.out, phenotype)
    }  
    LMM.confirmed %>% write_delim(res.file, delim=" ")
    cat(sprintf("List of verified LMM clumps written on:\n %s\n", res.file))
} else {
    # Load clumped discoveries
    clumps.file <- sprintf("%s/discoveries/%s_lmm_regions.txt", scratch, phenotype)
    LMM <- read_delim(clumps.file, delim=" ", col_types=cols())
    # Try to replicate
    LMM.confirmed <- LMM %>% verify_with_meta(Meta.pvalues, significance=parse_number(significance.meta))
    cat(sprintf("Verified LMM discoveries: %d out of %d (%.2f%%).\n",
            sum(!is.na(LMM.confirmed$Verified.BP)),
            nrow(LMM.confirmed),
            100*mean(!is.na(LMM.confirmed$Verified.BP))))
    # Save results
    res.file <- sprintf("%s/lmm_%s.txt", scratch.out, phenotype)
    LMM.confirmed %>% write_delim(res.file, delim=" ")
    cat(sprintf("List of verified LMM clumps written on:\n %s\n", res.file))
}

###################################
## Load clumped meta discoveries ##
###################################

# Load raw clumped results
clumped.file <- sprintf("%s/meta/clumped/%s_%s.clumped", scratch, phenotype.code, significance.meta)
Meta.raw <- read.table(clumped.file, header=TRUE) %>% as_tibble() %>%
    select(CHR, SNP, BP, SP2)
cat(sprintf("Loaded %d clumps for %s.\n", nrow(Meta.raw), phenotype))

# Parse clumped results
if(nrow(Meta.raw)>0) {
    # Extract lead SNPs
    Meta.raw.lead <- Meta.raw %>% mutate(SP2=SNP) %>% select(CHR, SNP, BP, SP2)
    # Extract secondary SNPs
    Meta.raw.secondary <- Meta.raw %>%
        mutate(SP2 = gsub("\\(1\\)", "", SP2)) %>%
        separate_rows(SP2, sep=",") %>%
        select(CHR, SNP, BP, SP2) %>%
        filter(SP2!="NONE")
    # Cross-reference with complete list of variants
    Meta <- rbind(Meta.raw.lead, Meta.raw.secondary) %>%
        as_tibble() %>% arrange(CHR, BP, SP2) %>%
        mutate(SNP.lead=as.character(SNP), BP.lead=BP, SNP=as.character(SP2)) %>%
        select(CHR, SNP.lead, BP.lead, SNP) %>%
        left_join(Meta.pvalues, by=c('CHR', 'SNP')) %>%
        select(CHR, SNP.lead, BP.lead, SNP, BP, P)
} else {
    Meta <- tibble()
}

Meta <- Meta %>%
    group_by(CHR, SNP.lead, BP.lead) %>%
    summarise(BP.min=min(BP), BP.max=max(BP), P=min(P)) %>%
    ungroup()

# Consolidate the discoveries
Meta.consolidated <- Meta %>% mutate(Size=NA) %>% consolidate_clumps()

####################################
## Meta discoveries vs. knockoffs ##
####################################
match_to_discoveries <- function(chr, bp.min, bp.max, Discoveries) {
    if("P" %in% colnames(Discoveries)) {
        Discoveries <- Discoveries %>% mutate(Importance=-log10(P))
    } else {
        Discoveries <- Discoveries %>% mutate(Importance=W)
    }
    gap <- 1e5
    matches <- Discoveries %>%
        filter(CHR==chr, BP.min<=bp.max+gap, BP.max>=bp.min-gap) %>%
        arrange(desc(Importance))
    if(nrow(matches)>0) {
        return(matches$SNP.lead[1])
    } else {
        return(NA)
    }
    return(matches$BP[1])
}
verify_with_clumps <- function(Clumped, Discoveries) {
    # Match each row of the discoveries with the meta analysis
    Clumped <- Clumped %>%
        rowwise() %>%
        mutate(Verified.SNP = match_to_discoveries(CHR, BP.min, BP.max, Discoveries)) %>%
        ungroup()
    # Add BP info
    Discoveries.ref <- Discoveries %>%
        mutate(Verified.BP=BP.lead, Verified.SNP=SNP.lead) %>%
        select(CHR, Verified.BP, Verified.SNP)
    Clumped <- Clumped %>%
        left_join(Discoveries.ref, by=c("CHR", "Verified.SNP"))
    return(Clumped)
}

# Confirm discoveries (consolidated)
Meta.confirmed <- Meta.consolidated %>% verify_with_clumps(Knockoffs)
cat(sprintf("Consolidated Meta discoveries verified by knockoffs: %d out of %d (%.2f%%).\n",
            sum(!is.na(Meta.confirmed$Verified.SNP)),
            nrow(Meta.consolidated),
            100*mean(!is.na(Meta.confirmed$Verified.SNP))
            ))
# Save results
if(use.small) {
    res.file <- sprintf("%s/meta_consolidated_knockoffs_%s_%s_fold_01.txt", scratch.out, phenotype, significance.meta)
} else {
    res.file <- sprintf("%s/meta_consolidated_knockoffs_%s_%s.txt", scratch.out, phenotype, significance.meta)
}
Meta.confirmed %>% write_delim(res.file, delim=" ")
cat(sprintf("List of verified LMM clumps (full data, consolidated) written on:\n %s\n", res.file))

# Confirm discoveries (not consolidated)
Meta.confirmed <- Meta %>% verify_with_clumps(Knockoffs)
cat(sprintf("Consolidated Meta discoveries verified by knockoffs: %d out of %d (%.2f%%).\n",
            sum(!is.na(Meta.confirmed$Verified.SNP)),
            nrow(Meta),
            100*mean(!is.na(Meta.confirmed$Verified.SNP))
            ))
# Save results
if(use.small) {
    res.file <- sprintf("%s/meta_knockoffs_%s_%s_fold_01.txt", scratch.out, phenotype, significance.meta)
} else {
    res.file <- sprintf("%s/meta_knockoffs_%s_%s.txt", scratch.out, phenotype, significance.meta)
}
Meta.confirmed %>% write_delim(res.file, delim=" ")
cat(sprintf("List of verified LMM clumps (full data) written on:\n %s\n", res.file))

####################################
## Meta discoveries vs. LMM       ##
####################################

# Confirm discoveries (consolidated)
Meta.confirmed <- Meta.consolidated %>% verify_with_clumps(LMM)
cat(sprintf("Meta discoveries verified by LMM: %d out of %d (%.2f%%).\n",
            sum(!is.na(Meta.confirmed$Verified.SNP)),
            nrow(Meta.consolidated),
            100*mean(!is.na(Meta.confirmed$Verified.SNP))
            ))
# Save results
if(use.small) {
    if(use.BH) {        
        res.file <- sprintf("%s/meta_consolidated_lmm_%s_%s_fold_01_BH.txt",
                            scratch.out, phenotype, significance.meta)
    } else {
        res.file <- sprintf("%s/meta_consolidated_lmm_%s_%s_fold_01.txt",
                            scratch.out, phenotype, significance.meta)
    }
} else {
    res.file <- sprintf("%s/meta_consolidated_lmm_%s_%s.txt", scratch.out, phenotype, significance.meta)
}
Meta.confirmed %>% write_delim(res.file, delim=" ")
cat(sprintf("List of verified LMM clumps (full data, consolidated) written on:\n %s\n", res.file))

# Confirm discoveries (not consolidated)
Meta.confirmed <- Meta %>% verify_with_clumps(LMM)
cat(sprintf("Meta discoveries verified by LMM: %d out of %d (%.2f%%).\n",
            sum(!is.na(Meta.confirmed$Verified.SNP)),
            nrow(Meta),
            100*mean(!is.na(Meta.confirmed$Verified.SNP))
            ))
# Save results
if(use.small) {
    if(use.BH) {        
        res.file <- sprintf("%s/meta_lmm_%s_%s_fold_01_BH.txt", scratch.out, phenotype, significance.meta)
    } else {
        res.file <- sprintf("%s/meta_lmm_%s_%s_fold_01.txt", scratch.out, phenotype, significance.meta)
    }
} else {
    res.file <- sprintf("%s/meta_lmm_%s_%s.txt", scratch.out, phenotype, significance.meta)
}
Meta.confirmed %>% write_delim(res.file, delim=" ")
cat(sprintf("List of verified LMM clumps (full data) written on:\n %s\n", res.file))

