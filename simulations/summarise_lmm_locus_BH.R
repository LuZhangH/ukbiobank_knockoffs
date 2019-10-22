#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

args <- commandArgs(trailingOnly=TRUE)

experiment.name <- "six"
experiment.name <- as.character(args[1])

library(tidyverse)
library(doParallel)
source("../utils/utils_clumping.R")

K <- 50
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1, 22)

amplitude.list <- seq(1,10)
fold.list <- seq(1,1)
resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")

###########################
## Load list of variants ##
###########################

resolution.dummy <- "Radj100"
Params <- expand.grid(Resolution=resolution.dummy, CHR=chr.list) %>%
    as_tibble() %>% mutate(Resolution=as.character(Resolution), CHR=as.integer(CHR))
Variants <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    resolution <- Params$Resolution[idx]
    chr <- Params$CHR[idx]
    cat(sprintf("Loading list of variants at resolution %s for chromosome %d...\n", resolution, chr))
    # Load list of variants
    key.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.key", scratch, resolution, K, chr)
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
    filter(Resolution==resolution.dummy, Knockoff==FALSE) %>%
    mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP))

###############################
## Load genetic architecture ##
###############################
arch.file <- sprintf("%s/simulations/phenotypes/architecture.txt", scratch)
Architecture <- read_delim(arch.file, delim=" ", col_types=cols()) %>%
    filter(Name==experiment.name)

##############################
## Load clumped discoveries ##
##############################

Params <- expand.grid(Amplitude=amplitude.list, Fold=fold.list) %>% as_tibble()
Params.loaded <- Params %>% head(0)

Selections <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    amplitude <- Params$Amplitude[idx]
    fold <- Params$Fold[idx]
    cat(sprintf("Loading discoveries for fold %s, amplitude %d...\n", fold, amplitude))
    # Load clumped discoveries
    clumped.file <- sprintf("%s/simulations/BH/fold_1/selections_%s_%s.txt",
                            scratch, experiment.name, amplitude)
    if(file.exists(clumped.file)) {        
        Params.loaded <- rbind(Params.loaded, Params[idx,])
        Clumped <- read_delim(clumped.file, delim=" ", col_types=cols())       
        # Add meta-information and return list of discoveries
        if(nrow(Clumped)>0) {
            Clumped <- cbind(Params[idx,], Clumped) %>% as_tibble()
        }
    } else {
        Clumped <- tibble()
    }
    return(Clumped)
}

loci.intersect <- function(bp.min.x, bp.max.x, bp.min.y, bp.max.y, gap=1e5) {
    intersect <- as.logical((bp.min.y<=bp.max.x+gap) * (bp.max.y>=bp.min.x-gap))
    return(intersect)
}

##################################
## Summarise results (distinct) ##
##################################

Selections.loci <- Selections %>%
    full_join(Params, by = colnames(Params)) %>%
    full_join(Architecture %>% select(Amplitude, H2), by = "Amplitude") %>%
    group_by(Amplitude, H2, Fold, Resolution, CHR, Group) %>%
    summarise(P=min(P.simes), SNP.lead=SNP[1], BP.lead=BP[1], BP.min=min(BP), BP.max=max(BP),
              Size=n(), Causal=any(Causal)) %>%
    mutate(BP.width=BP.max-BP.min) %>%
    ungroup() %>%
    mutate(Resolution=factor(Resolution, levels=resolution.list, labels=resolution.list))

Distinct.summary <- Selections.loci %>% group_by(Amplitude, H2, Fold, Resolution) %>%
    summarise(Discoveries=sum(!is.na(Group)), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              Size=mean(Size),BP.width=mean(BP.width)) %>%
    mutate(Method="LMM-BH")

cat("Summary of distinct discoveries with LMM + BH:\n")
Distinct.summary %>% print()

out.file <- sprintf("%s/simulations/summary/lmm_BH_%s_distinct.txt", scratch, experiment.name)
Distinct.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

Selections.loci.consolidated <- Selections.loci %>%
    filter(Resolution=="Radj2") %>%
    group_by(Amplitude, H2, Fold, Resolution) %>%
    consolidate_clumps(gap=1e5)

Distinct.consolidated.summary <- Selections.loci.consolidated %>%
    group_by(Amplitude, H2, Fold, Resolution) %>%
    summarise(Discoveries=sum(!is.na(SNP.lead)), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              BP.width=mean(BP.width), Size=mean(Size)) %>%
    mutate(Method="LMM-BH")

cat("Summary of distinct consolidated discoveries with LMM + BH:\n")
Distinct.consolidated.summary %>% print()

out.file <- sprintf("%s/simulations/summary/lmm_BH_%s_distinct_consolidated.txt", scratch, experiment.name)
Distinct.consolidated.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

#####################################
## Summarise results (locus power) ##
#####################################

# Define causal clusters
Causal.clusters <- Causal %>%
    mutate(BP.lead=BP, BP.min=BP,BP.max=BP,Importance=1,SNP.lead=BP,Size=1) %>%
    consolidate_clumps(gap=1e5) %>%
    select(CHR, BP.lead, BP.max, BP.min, Size)

# Evaluate power to recover clusters
Cluster.summary <- foreach(idx = 1:nrow(Params.loaded), .combine = 'rbind') %dopar% {
    amplitude <- Params.loaded$Amplitude[idx]
    fold <- Params.loaded$Fold[idx]
    Loci <- Selections.loci.consolidated %>%
        filter(Resolution=="Radj2", Amplitude==amplitude, Fold==fold)
    cat(sprintf("Computing power at the locus level for amplitude %d and fold %s...\n", amplitude, fold))
    # Compute power
    if(nrow(Loci)>0) {
        Causal.pow <- Causal.clusters %>% select(CHR, BP.min, BP.max) %>%
            left_join(Loci %>% select(CHR, BP.min, BP.max), by="CHR") %>%
            mutate(Intersect = loci.intersect(BP.min.x, BP.max.x, BP.min.y, BP.max.y)) %>%
            group_by(CHR, BP.min.x, BP.max.x) %>%
            summarise(Discovered=any(Intersect)) %>%
            ungroup() %>%
            mutate(Discovered = replace_na(Discovered, FALSE))
        pow <- mean(Causal.pow$Discovered)
    } else {
        pow <- 0
    }
    # Compute FDP
    if(nrow(Loci)>0) {
        Causal.fdp <- Loci %>% select(CHR, BP.min, BP.max) %>%
            left_join(Causal.clusters %>% select(CHR, BP.min, BP.max), by="CHR") %>%
            mutate(Intersect = loci.intersect(BP.min.x, BP.max.x, BP.min.y, BP.max.y)) %>%
            group_by(CHR, BP.min.x, BP.max.x) %>%
            summarise(Causal=any(Intersect))
        fdp <- mean(!Causal.fdp$Causal)
    } else {
        fdp <- 0
    }
    Results.row <- tibble(Amplitude=amplitude, Fold=fold, Method="LMM-BH",
                          Discoveries=sum(!is.na(Loci$SNP.lead)),
                          Discoveries.true=sum(Loci$Causal,na.rm=T), Power=pow, FDP=fdp,
                          Size=mean(Loci$Size), BP.width=mean(Loci$BP.width))
}

Cluster.summary <- Cluster.summary %>% inner_join(Architecture %>% select(Amplitude, H2), by=c("Amplitude"))

cat("Summary of locus discovery with LMM:\n")
Cluster.summary %>% print()

out.file <- sprintf("%s/simulations/summary/lmm_BH_%s_loci.txt", scratch, experiment.name)
Cluster.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))
