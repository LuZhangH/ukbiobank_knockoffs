#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")
library(tidyverse)
source("../utils/utils_clumping.R")

phenotype.code <- "body_HEIGHTz" # "bp_SYSTOLICadjMEDz" "body_HEIGHTz"
    
args <- commandArgs(trailingOnly=TRUE)
phenotype.code <- as.character(args[1])

significance.meta <- "0.00000005" # 5e-8
resolution <- "Radj2"
offset <- 1

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
Knockoffs <- lapply(c(0,1), function(seed) {
  discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_%s_s%s_offset%d.txt",
                              scratch, phenotype, resolution, seed, offset)
  Knockoffs.s <- read_delim(discoveries.file, delim=" ", col_types=cols()) %>% mutate(Seed=seed)
  return(Knockoffs.s)
})
Knockoffs <- do.call("rbind", Knockoffs)

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

cat(sprintf("Summary of knockoff discoveries:\n"))
Knockoffs.confirmed %>%
  group_by(Seed) %>% summarise(Discoveries=n(), Confirmed=sum(!is.na(Verified.SNP))) %>%
  print()

## Compare discoveries obtained with different seeds
Knockoffs.0 <- Knockoffs.confirmed %>% filter(Seed==0) %>% select(-Seed)
Knockoffs.1 <- Knockoffs.confirmed %>% filter(Seed==1) %>% select(-Seed)
vars.joinby <- setdiff(colnames(Knockoffs.0), "W")
Discoveries.all <- full_join(Knockoffs.0, Knockoffs.1, by=vars.joinby)
Knockoffs.both <- Discoveries.all %>% filter(!is.na(W.x), !is.na(W.y))

Discoveries.all %>%
    summarise(N0 = sum(!is.na(W.x)), N1 = sum(!is.na(W.y)),
              N = sum((!is.na(W.x))*(!is.na(W.y)))) %>%
    mutate(Stable=N/pmax(1,N0)*100)

Discoveries.all %>% mutate(Confirmed=!is.na(Verified.SNP)) %>%
  group_by(Confirmed) %>%
  summarise(N0 = sum(!is.na(W.x)), N1 = sum(!is.na(W.y)),
            N = sum((!is.na(W.x))*(!is.na(W.y)))) %>%
  mutate(Stable=N/pmax(1,N0)*100)

res.file <- sprintf("%s/meta_knockoffs_stability_%s.txt", scratch.out, phenotype)
Knockoffs.confirmed %>% write_delim(res.file, delim=" ")
cat(sprintf("List of verified knockoff discoveries written on:\n %s\n", res.file))

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

#######################
## Meta VS Knockoffs ##
#######################

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
Meta.confirmed.0 <- Meta.consolidated %>% verify_with_clumps(Knockoffs.0) %>% mutate(Seed=0)
Meta.confirmed.1 <- Meta.consolidated %>% verify_with_clumps(Knockoffs.1) %>% mutate(Seed=1)  
Meta.confirmed.both <- inner_join(Meta.confirmed.0, Meta.confirmed.1,
                                  by = c("CHR", "SNP.lead", "BP.lead", "BP.min", "BP.max", "P", "Size")) %>%
  mutate(Verified.SNP=ifelse((!is.na(Verified.SNP.x))*(!is.na(Verified.SNP.x)), Verified.SNP.x, NA)) %>%
  mutate(Verified.BP=ifelse((!is.na(Verified.BP.x))*(!is.na(Verified.BP.x)), Verified.BP.x, NA)) %>%
  select(-Verified.SNP.x, -Verified.SNP.y, -Verified.BP.x, -Verified.BP.y, -Seed.x, -Seed.y) %>%
  mutate(Seed="both")

Meta.confirmed <- rbind(Meta.confirmed.0, Meta.confirmed.1, Meta.confirmed.both)

res.file <- sprintf("%s/meta_lmm_consolidated_stability_%s.txt", scratch.out, phenotype)
Meta.confirmed %>% write_delim(res.file, delim=" ")
cat(sprintf("List of verified LMM discoveries written on:\n %s\n", res.file))

cat(sprintf("Consolidated Meta discoveries verified by knockoffs (seed %s): %d out of %d (%.2f%%).\n",
            0,
            sum(!is.na(Meta.confirmed.0$Verified.SNP)),
            nrow(Meta.consolidated),
            100*mean(!is.na(Meta.confirmed.0$Verified.SNP))
            ))
cat(sprintf("Consolidated Meta discoveries verified by knockoffs (seed %s): %d out of %d (%.2f%%).\n",
            1,
            sum(!is.na(Meta.confirmed.1$Verified.SNP)),
            nrow(Meta.consolidated),
            100*mean(!is.na(Meta.confirmed.1$Verified.SNP))
            ))
cat(sprintf("Consolidated Meta discoveries verified by knockoffs (seed %s): %d out of %d (%.2f%%).\n",
            "both",
            sum(!is.na(Meta.confirmed.both$Verified.SNP)),
            nrow(Meta.consolidated),
            100*mean(!is.na(Meta.confirmed.both$Verified.SNP))
            ))

# Confirm discoveries (not consolidated)
Meta.confirmed.0 <- Meta %>% verify_with_clumps(Knockoffs.0) %>% mutate(Seed=0)
Meta.confirmed.1 <- Meta %>% verify_with_clumps(Knockoffs.1) %>% mutate(Seed=1)
Meta.confirmed.both <- inner_join(Meta.confirmed.0, Meta.confirmed.1,
                                  by = c("CHR", "SNP.lead", "BP.lead", "BP.min", "BP.max", "P")) %>%
  mutate(Verified.SNP=ifelse((!is.na(Verified.SNP.x))*(!is.na(Verified.SNP.x)), Verified.SNP.x, NA)) %>%
  mutate(Verified.BP=ifelse((!is.na(Verified.BP.x))*(!is.na(Verified.BP.x)), Verified.BP.x, NA)) %>%
  select(-Verified.SNP.x, -Verified.SNP.y, -Verified.BP.x, -Verified.BP.y, -Seed.x, -Seed.y) %>%
  mutate(Seed="both")

Meta.confirmed <- rbind(Meta.confirmed.0, Meta.confirmed.1, Meta.confirmed.both)

res.file <- sprintf("%s/meta_lmm_stability_%s.txt", scratch.out, phenotype)
Meta.confirmed %>% write_delim(res.file, delim=" ")
cat(sprintf("List of verified LMM discoveries written on:\n %s\n", res.file))

cat(sprintf("Meta discoveries verified by knockoffs (seed %s): %d out of %d (%.2f%%).\n",
            0,
            sum(!is.na(Meta.confirmed.0$Verified.SNP)),
            nrow(Meta),
            100*mean(!is.na(Meta.confirmed.0$Verified.SNP))
            ))
cat(sprintf("Meta discoveries verified by knockoffs (seed %s): %d out of %d (%.2f%%).\n",
            1,
            sum(!is.na(Meta.confirmed.1$Verified.SNP)),
            nrow(Meta),
            100*mean(!is.na(Meta.confirmed.1$Verified.SNP))
            ))
cat(sprintf("Meta discoveries verified by knockoffs (seed %s): %d out of %d (%.2f%%).\n",
            "both",
            sum(!is.na(Meta.confirmed.both$Verified.SNP)),
            nrow(Meta),
            100*mean(!is.na(Meta.confirmed.both$Verified.SNP))
            ))

