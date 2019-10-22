#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

experiment.name <- "six"
conservative <- 0

args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])
conservative <- as.integer(args[2])

library(tidyverse)
library(purrrlyr)
library(doParallel)
source("../utils/utils_clumping.R")
source("utils_outer.R")

K <- 50
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1, 22)
resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")
#resolution.list <- c("Radj2")

amplitude.list <- seq(1,10)
fold.list <- seq(1,1)

###############################
## Load genetic architecture ##
###############################
arch.file <- sprintf("%s/simulations/phenotypes/architecture.txt", scratch)
Architecture <- read_delim(arch.file, delim=" ", col_types=cols()) %>%
    filter(Name==experiment.name)

# Load list of causal variants
causal.file <- sprintf("%s/simulations/phenotypes/%s_causal.txt", scratch, experiment.name)
Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, BP, Locus, Sign, Scale)

##############################
## Load knockoff statistics ##
##############################

Params <- expand.grid(Amplitude=amplitude.list, Fold=fold.list) %>% as_tibble()

Selections <- lapply(1:nrow(Params), function(idx) {
    fold <- Params$Fold[idx]
    amplitude <- Params$Amplitude[idx]
    cat(sprintf("Loading knockoff statistics for fold %s, amplitude %d...\n", fold, amplitude))
    ## Load discoveries
    discoveries.file <- sprintf("%s/simulations/lasso/fold_%s/multires/lasso_%s_%s_c%s.txt",
                            scratch, fold, experiment.name, amplitude, conservative)
    if(file.exists(discoveries.file)) {
      Discoveries <- read_delim(discoveries.file, delim=" ", col_types=cols()) %>%
        mutate(Fold=fold, Amplitude=amplitude)
    } else {
      warning(sprintf("File %s does not exist.\n", discoveries.file))
      Discoveries <- tibble()
    }
    return(Discoveries)
})
Selections <- do.call("rbind", Selections)

# Preview of solutions
cat(sprintf("Solutions at resolution Radj2:\n"))
Selections %>% filter(Resolution=="Radj2", Fold==1) %>%
    group_by(Amplitude, Resolution, Fold) %>%
    summarise(N=n(), True=sum(Causal), False=sum(!Causal), FDP=mean(!Causal)) %>%
    ungroup() %>%
    group_by(Amplitude, Resolution) %>%
    summarise(FDP=mean(FDP))

##########################
## Outer-node knockoffs ##
##########################

Params.outer <- expand.grid(Fold=fold.list, Amplitude=amplitude.list) %>%
    as_tibble() %>% mutate(Fold=as.character(Fold))

Selections.outer <- lapply(1:nrow(Params.outer), function(idx) {
    fold <- Params$Fold[idx]
    amplitude <- Params$Amplitude[idx]
    Selections.outer.row <- Selections %>% mutate(Resolution=parse_number(Resolution)) %>%
        filter(Amplitude==amplitude, Fold==fold)
    if(nrow(Selections.outer.row)>0) {
        cat(sprintf("Resolution filter for amplitude %d and fold %s ... \n", amplitude, fold))
        Selections.outer.row <- Selections.outer.row %>%
            consistent_filter() %>% consistent_to_outer()
        Selections.outer.row$Resolution <- "outer"
    }
    return(Selections.outer.row)
})
Selections.outer <- do.call("rbind", Selections.outer)

cat(sprintf("Outer node selections:\n"))
Selections.outer %>% group_by(Amplitude) %>%
    summarise(N=n(), True=sum(Causal), False=sum(!Causal), FDP=mean(!Causal)) %>% print()
Selections <- rbind(Selections, Selections.outer)

Params <- expand.grid(Amplitude=amplitude.list, Fold=fold.list,
                      Resolution=c(resolution.list,"outer")) %>%
    as_tibble() %>%
    mutate(Resolution=as.character(Resolution))

##################################
## Summarise results (distinct) ##
##################################

Selections.loci <- Selections %>%
    full_join(Params, by = colnames(Params)) %>%
    full_join(Architecture %>% select(Amplitude, H2), by = "Amplitude") %>%
    ungroup()

Distinct.summary <- Selections.loci %>% group_by(Amplitude, H2, Fold, Resolution) %>%
    summarise(Discoveries=n(), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              BP.width=mean(BP.width), Size=mean(Size)) %>%
    mutate(Method="Knockoffs")

cat("Summary of distinct discoveries with Knockoffs:\n")
Distinct.summary %>% group_by(Amplitude, H2, Resolution) %>%
    summarise(Discoveries.true=mean(Discoveries.true), FDP=mean(FDP)) %>%
    filter(Resolution=="Radj2") %>% print()

out.file <- sprintf("%s/simulations/summary/knockoffs_%s_distinct_multires_c%s.txt",
                    scratch, experiment.name, conservative)
Distinct.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

Selections.loci.consolidated <- Selections.loci %>%
    mutate(SNP.lead=NA, BP.lead=NA, Importance=1) %>%
    group_by(Amplitude, H2, Fold, Resolution) %>%
    consolidate_clumps(gap=1e5)

Distinct.consolidated.summary <- Selections.loci.consolidated %>%
    group_by(Amplitude, H2, Fold, Resolution) %>%
    summarise(Discoveries=sum(!is.na(SNP.lead)), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              BP.width=mean(BP.width), Size=mean(Size)) %>%
    mutate(Method="Knockoffs")

cat("Summary of distinct consolidated discoveries with knockoffs:\n")
Distinct.consolidated.summary %>% print()

out.file <- sprintf("%s/simulations/summary/knockoffs_%s_distinct_consolidated_multires_c%s.txt",
                    scratch, experiment.name, conservative)
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
Cluster.summary <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    amplitude <- Params$Amplitude[idx]
    fold <- Params$Fold[idx]
    resolution <- Params$Resolution[idx]
    Loci <- Selections.loci.consolidated %>%
        filter(Amplitude==amplitude, Resolution==resolution, Fold==fold)
    cat(sprintf("Computing power at the locus level for amplitude %d and resolution %s, fold %s...\n",
                amplitude, resolution, fold))
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
    Results.row <- tibble(Amplitude=amplitude, Fold=fold, Method="Knockoffs", Resolution=resolution,
                          Discoveries=sum(!is.na(Loci$SNP.lead)), Discoveries.true=sum(Loci$Causal,na.rm=T),
                          Power=pow, FDP=fdp,
                          Size=mean(Loci$Size), BP.width=mean(Loci$BP.width))
}

Cluster.summary <- Cluster.summary %>%
    inner_join(Architecture %>% select(Amplitude, H2), by=c("Amplitude"))

cat("Summary of locus discovery with Knockoffs:\n")
Cluster.summary %>% print()

out.file <- sprintf("%s/simulations/summary/knockoffs_%s_loci_multires_c%s.txt",
                    scratch, experiment.name, conservative)
Cluster.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))
