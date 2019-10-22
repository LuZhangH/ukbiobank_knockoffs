#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

# Main parameters (default)
experiment.name <- "seven"
amplitude <- 10

# Main parameters (input)
args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])
amplitude <- as.numeric(args[2])

library(tidyverse)
library(doParallel)

K <- 50
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1, 22)

###############################
## Load genetic architecture ##
###############################

arch.file <- sprintf("%s/simulations/phenotypes/architecture.txt", scratch)
Architecture <- read_delim(arch.file, delim=" ", col_types=cols()) %>%
    filter(Name==experiment.name)

####################################
## Load list of original variants ##
####################################
resolution.dummy <- "Radj100"
Params <- expand.grid(CHR=chr.list) %>%
    as_tibble() %>% mutate(CHR=as.integer(CHR))
Variants <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    chr <- Params$CHR[idx]
    cat(sprintf("Loading list of variants at resolution %s for chromosome %d...\n", resolution.dummy, chr))
    # Load list of variants
    key.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.key", scratch, resolution.dummy, K, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants.chr$Resolution <- resolution.dummy
    return(Variants.chr)
}
# Load list of causal variants
causal.file <- sprintf("%s/simulations/phenotypes/%s_causal.txt", scratch, experiment.name)
Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, BP, Locus, Sign, Scale)
Variants <- Variants %>%
    left_join(mutate(Causal,Causal=TRUE),by = c("CHR", "BP")) %>%
    mutate(Causal=replace_na(Causal,FALSE)) %>%
    mutate(Causal=ifelse(Knockoff, FALSE, Causal))
# Extract list of original variants
Variants.original <- Variants %>%
    filter(Knockoff==FALSE) %>%
    mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP))

###################
## Load p-values ##
###################

lmm.file <- sprintf("%s/simulations/bolt/fold_1/stats_%s_%s.txt", scratch, experiment.name, amplitude)
LMM.raw <- read_delim(lmm.file, delim="\t", col_types=cols())
if("P_BOLT_LMM" %in% colnames(LMM.raw)) {
    LMM.raw <- LMM.raw %>% mutate(P=P_BOLT_LMM)
} else {
    LMM.raw <- LMM.raw %>% mutate(P=P_BOLT_LMM_INF)
}

#############################################
## Load list of variants and select groups ##
#############################################

resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")

Results <- lapply(resolution.list, function(resolution) {

    Params <- expand.grid(CHR=chr.list) %>%
        as_tibble() %>% mutate(CHR=as.integer(CHR))
    Variants <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
        chr <- Params$CHR[idx]
        cat(sprintf("Loading list of variants at resolution %s for chromosome %d...\n", resolution, chr))
        # Load list of variants
        key.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.key", scratch, resolution, K, chr)
        Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
        Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
        colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
        Variants.chr$Resolution <- resolution
        return(Variants.chr)
    }
    # Load list of causal variants
    causal.file <- sprintf("%s/simulations/phenotypes/%s_causal.txt", scratch, experiment.name)
    Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, BP, Locus, Sign, Scale)
    Variants <- Variants %>%
        left_join(mutate(Causal,Causal=TRUE),by = c("CHR", "BP")) %>%
        mutate(Causal=replace_na(Causal,FALSE)) %>%
        mutate(Causal=ifelse(Knockoff, FALSE, Causal))
    # Extract list of original variants
    Variants.original <- Variants %>%
        filter(Resolution==resolution, Knockoff==FALSE) %>%
        mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP))

    # Add variant information to the LMM p-values
    LMM <- LMM.raw %>% inner_join(Variants.original, by = c("SNP", "CHR", "BP"))
    
    # Summarise p-values by group
    LMM.groups <- LMM %>%
        group_by(CHR, Group) %>%
        summarise(N=n(), P.simes=min(sort(P)*N/seq(N)), P.bonferroni=min(1,min(P)*N)) %>%
        select(CHR, Group, N, P.simes, P.bonferroni)

    # Apply BH
    Selected <- LMM.groups %>%
        ungroup() %>%
        mutate(P.adjusted=p.adjust(P.simes, method="BH")) %>%
        filter(P.adjusted<=0.1)

    # Count discoveries
    Results.res <- Selected %>%
        inner_join(Variants.original, by=c("CHR", "Group")) %>%
        mutate(Resolution=resolution)

    # Give preview of results
    cat(sprintf("Results for resolution %s:\n", resolution))
    Results.res %>%
        group_by(CHR, Group) %>%
        summarise(Causal=any(Causal)) %>%
        ungroup() %>%
        summarise(Discoveries=n(), True=sum(Causal), FDP=mean(!Causal)) %>%
        print()
    
    return(Results.res)
})
Results <- do.call("rbind", Results)

##########################
## Evaluate performance ##
##########################

# Give preview of results
cat(sprintf("Results for all resolutions:\n"))
Summary <- Results %>%
    group_by(Resolution, CHR, Group) %>%
    summarise(Causal=any(Causal)) %>%
    ungroup() %>% group_by(Resolution) %>%
    summarise(Discoveries=n(), True=sum(Causal), FDP=mean(!Causal))
Summary %>% print()

# Save results
out.file <- sprintf("%s/simulations/BH/fold_1/selections_%s_%s.txt", scratch, experiment.name, amplitude)
Results %>%
    mutate(Resolution=factor(Resolution, levels=resolution.list, labels=resolution.list)) %>%
    arrange(Resolution, CHR, Group) %>%
    write_delim(out.file, delim=" ")
cat(sprintf("Results saved on:\n %s\n", out.file))
