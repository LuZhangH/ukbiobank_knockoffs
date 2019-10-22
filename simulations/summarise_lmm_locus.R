#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

experiment.name <- "six"

args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])

library(tidyverse)
library(doParallel)
source("../utils/utils_clumping.R")

K <- 50
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1, 22)
resolution.dummy <- "Radj50"

amplitude.list <- seq(1,10)
#amplitude.list <- seq(8,8)
fold.list <- seq(1,1)
clumping.list <- c("0.000000005","0.00000005","0.000005","0.00001","0.00002","0.00005","0.0001","0.0002","0.0005","0.001", "BH")
#clumping.list <- c("0.00000005")
#clumping.list <- c("0.00000005", "BH")

###########################
## Load list of variants ##
###########################

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

Params <- expand.grid(Amplitude=amplitude.list, Fold=fold.list, Clumping=clumping.list) %>% as_tibble() %>%
    mutate(Clumping=as.character(Clumping))
Params.loaded <- Params %>% head(0)

Selections <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    amplitude <- Params$Amplitude[idx]
    fold <- Params$Fold[idx]
    clumping <- Params$Clumping[idx]
    cat(sprintf("Loading discoveries for fold %s, amplitude %d, clumping %s...\n",
                fold, amplitude, clumping))

    # Load clumped discoveries
    clumped.file <- sprintf("%s/simulations/bolt_clumped/fold_%s/stats_%s_%s_%s.clumped",
                            scratch, fold, experiment.name, amplitude, clumping)
    if(file.exists(clumped.file)) {
        Params.loaded <- rbind(Params.loaded, Params[idx,])
        Clumped.raw <- read.table(clumped.file, header=TRUE) %>% as_tibble() %>%
            select(CHR, SNP, BP, P, SP2)

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

########################
## Oracle calibration ##
########################

calibrate_oracle <- function(fdr, power, fdr.target=0.1) {
    if(length(fdr)==1) {
        return(power)
    }
    power <- c(0,power)
    fdr <- c(0,fdr)
    fdr[is.na(fdr)] <- 0

    func <- splinefun(x=fdr, y=power, method="monoH.FC",  ties = mean)
    power.interpolate <- func(0.1)

    return(power.interpolate)
}

##################################
## Summarise results (distinct) ##
##################################

Selections.loci <- Selections %>%
    full_join(Params, by = colnames(Params)) %>%
    full_join(Architecture %>% select(Amplitude, H2), by = "Amplitude") %>%
    group_by(Amplitude, H2, Fold, Clumping, CHR, SNP.lead, BP.lead) %>%
    summarise(P=min(P), BP.min=min(BP), BP.max=max(BP), Size=n(), Causal=any(Causal)) %>%
    mutate(BP.width=BP.max-BP.min) %>%
    ungroup()

Distinct.summary <- Selections.loci %>% group_by(Amplitude, H2, Fold, Clumping) %>%
    summarise(Discoveries=sum(!is.na(SNP.lead)), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              Size=mean(Size),BP.width=mean(BP.width)) %>%
    mutate(Method="LMM")

Oracle <- Distinct.summary %>%
    filter(Clumping!="BH") %>%
    group_by(Method, Amplitude, H2, Clumping) %>%
    summarise(FDP=mean(FDP), BP.width=mean(BP.width), Size=mean(Size),
              Discoveries=mean(Discoveries), Discoveries.true=mean(Discoveries.true)) %>%
    group_by(Method, Amplitude, H2) %>%
    summarise(Discoveries=calibrate_oracle(FDP, Discoveries),
              Discoveries.true=calibrate_oracle(FDP, Discoveries.true),
              BP.width=calibrate_oracle(FDP, BP.width), Size=calibrate_oracle(FDP, Size),
              FDP=0.1, Clumping="oracle", Fold=0)
Distinct.summary <- rbind(Distinct.summary, Oracle)

cat("Summary of distinct discoveries with LMM:\n")
Distinct.summary %>% print()

out.file <- sprintf("%s/simulations/summary/lmm_%s_distinct.txt", scratch, experiment.name)
Distinct.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

Selections.loci.consolidated <- Selections.loci %>%
    group_by(Amplitude, H2, Fold, Clumping) %>%
    consolidate_clumps(gap=1e5)

Distinct.consolidated.summary <- Selections.loci.consolidated %>%
    group_by(Amplitude, H2, Fold, Clumping) %>%
    summarise(Discoveries=sum(!is.na(SNP.lead)), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              BP.width=mean(BP.width), Size=mean(Size)) %>%
    mutate(Method="LMM")

Oracle <- Distinct.consolidated.summary %>%
    filter(Clumping!="BH") %>%
    group_by(Method, Amplitude, H2, Clumping) %>%
    summarise(FDP=mean(FDP), BP.width=mean(BP.width), Size=mean(Size),
              Discoveries=mean(Discoveries), Discoveries.true=mean(Discoveries.true)) %>%
    group_by(Method, Amplitude, H2) %>%
    summarise(Discoveries=calibrate_oracle(FDP, Discoveries),
              Discoveries.true=calibrate_oracle(FDP, Discoveries.true),
              BP.width=calibrate_oracle(FDP, BP.width),
              FDP=0.1, Clumping="oracle", Fold=0)
Distinct.consolidated.summary <- rbind(Distinct.consolidated.summary, Oracle)

cat("Summary of distinct consolidated discoveries with LMM:\n")
Distinct.consolidated.summary %>% print()

out.file <- sprintf("%s/simulations/summary/lmm_%s_distinct_consolidated.txt", scratch, experiment.name)
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
    clumping <- Params.loaded$Clumping[idx]
    Loci <- Selections.loci.consolidated %>%
        filter(Amplitude==amplitude, Clumping==clumping, Fold==fold)
    cat(sprintf("Computing power at the locus level for amplitude %d and clumping %s, fold %s...\n",
                amplitude, clumping, fold))
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
    Results.row <- tibble(Amplitude=amplitude, Fold=fold, Method="LMM", Clumping=clumping,
                          Discoveries=sum(!is.na(Loci$SNP.lead)),
                          Discoveries.true=sum(Loci$Causal,na.rm=T), Power=pow, FDP=fdp,
                          Size=mean(Loci$Size), BP.width=mean(Loci$BP.width))
}

Cluster.summary <- Cluster.summary %>% inner_join(Architecture %>% select(Amplitude, H2), by=c("Amplitude"))

Oracle <- Cluster.summary %>%
    filter(Clumping!="BH") %>%
    group_by(Method, Amplitude, H2, Clumping) %>%
    summarise(FDP=mean(FDP), BP.width=mean(BP.width), Size=mean(Size),
              Power=mean(Power), Discoveries=mean(Discoveries), Discoveries.true=mean(Discoveries.true)) %>%
    group_by(Method, Amplitude, H2) %>%
    summarise(Power=calibrate_oracle(FDP, Power),
              Discoveries=calibrate_oracle(FDP, Discoveries),
              Discoveries.true=calibrate_oracle(FDP, Discoveries.true),
              Size=round(calibrate_oracle(FDP, Size)),
              BP.width=calibrate_oracle(FDP, BP.width),
              FDP=0.1, Clumping="oracle", Fold=0) %>%
    ungroup()
Cluster.summary <- rbind(Cluster.summary, Oracle)

cat("Summary of locus discovery with LMM:\n")
Cluster.summary %>% print()

out.file <- sprintf("%s/simulations/summary/lmm_%s_loci.txt", scratch, experiment.name)
Cluster.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

if(FALSE) {

    Cluster.summary %>%
        ggplot(aes(x=H2, y=Discoveries, color=Clumping)) + geom_point() + geom_line()

    Cluster.summary %>%
        ggplot(aes(x=H2, y=FDP, color=Clumping)) +
        geom_point() + geom_line() +
        ylim(0,1) +
        geom_hline(yintercept=0.1) +
        theme_bw()

}
