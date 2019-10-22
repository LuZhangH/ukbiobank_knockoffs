#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

experiment.name <- "seven"

args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])

library(tidyverse)
library(doParallel)
source("../utils/utils_clumping.R")

K <- 50
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1, 22)
resolution.dummy <- "Radj50"

#amplitude.list <- seq(1,10)
amplitude.list <- seq(1,10)
fold.list <- seq(1,1)
#clumping.list <- c("0.000000005","0.00000005","0.000005","0.00001","0.00002","0.00005","0.0001","0.0002","0.0005","0.001")
clumping.list <- c("0.00000005")
coverage <- 0.9
susie.nobias <- TRUE

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
## SUSIE                    ##
##############################

Params <- expand.grid(Amplitude=amplitude.list, Fold=fold.list, Clumping=clumping.list) %>%
    as_tibble() %>%
    mutate(Clumping=as.character(Clumping), Clumping.num=as.numeric(Clumping))
Params.loaded <- Params %>% head(0)

Selections.susie <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    amplitude <- Params$Amplitude[idx]
    fold <- Params$Fold[idx]
    clumping <- Params$Clumping[idx]
    cat(sprintf("Loading discoveries for fold %s, amplitude %d, clumping %.3g...\n",
                fold, amplitude, as.numeric(clumping)))

    # Load clumped discoveries
    if(susie.nobias) {
        fine.file <- sprintf("%s/simulations/susie/fold_%s/%s_%s_%s_%s_nobias.txt",
                             scratch, fold, experiment.name, amplitude, clumping, coverage)
    } else {
        fine.file <- sprintf("%s/simulations/susie/fold_%s/%s_%s_%s_%s.txt",
                             scratch, fold, experiment.name, amplitude, clumping, coverage)
    }
    if(file.exists(fine.file)) {
        Clumped.raw <- read_delim(fine.file, delim=" ", col_types=cols())
        if(nrow(Clumped.raw)>0) {
            Clumped <- Clumped.raw %>%
                mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP)) %>%
                mutate(SNP.lead = gsub(".A", "", SNP.lead), SNP.lead = gsub(".B", "", SNP.lead)) %>%
                select(CHR, SNP.lead, BP.lead, SNP) %>%
                inner_join(select(Variants.original,CHR,SNP,BP,Causal), by = c("CHR", "SNP"))
            # Add meta-information and return list of discoveries
            Clumped <- cbind(Params[idx,], Clumped) %>% as_tibble()
        } else {
            Clumped <- tibble()
        }
    } else {
        Clumped <- tibble()
    }
    return(Clumped)
}

Selections.susie.distinct <- Selections.susie %>%
    full_join(Params, by = colnames(Params)) %>%
    full_join(Architecture %>% select(Amplitude, H2), by = "Amplitude") %>%
    group_by(Amplitude, H2, Fold, Clumping, Clumping.num, SNP.lead, BP.lead) %>%
    summarise(BP.min=min(BP), BP.max=max(BP), Size=n(), Causal=any(Causal)) %>%
    mutate(BP.width=BP.max-BP.min) %>%
    ungroup()

Summary.susie <- Selections.susie.distinct %>%
    group_by(Amplitude, H2, Fold, Clumping, Clumping.num) %>%
    summarise(Discoveries=sum(!is.na(SNP.lead)), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              Size=mean(Size),BP.width=mean(BP.width)) %>%
    select(Amplitude, Discoveries, Discoveries.true, FDP, BP.width, Size, everything(), -Clumping.num) %>%
    ungroup() %>%
    mutate(Coverage=coverage)
    
cat("Summary of fine-mapping discoveries with SUSIE:\n")
Summary.susie %>% print()

if(susie.nobias) {
    out.file <- sprintf("%s/simulations/summary/lmm_%s_susie_nobias.txt", scratch, experiment.name)
} else { 
    out.file <- sprintf("%s/simulations/summary/lmm_%s_susie.txt", scratch, experiment.name)
}       
Summary.susie %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

##############################
## CAVIAR                   ##
##############################
load.caviar <- function(amplitude, fold, clumping, coverage=0.9) {
    posterior.file <- sprintf("%s/simulations/caviar/fold_%s/%s_%s_%s.txt",
                    scratch, fold, experiment.name, amplitude, clumping)
    Posterior <- tibble()
    if(file.exists(posterior.file)) {
        Posterior <- read_delim(posterior.file, delim=" ", col_types=cols())
    } else {
        Posterior <- tibble()
    }
    if(nrow(Posterior)>0) {
        # Select subset of variants based on CAVIAR posterior probabilities
        if(coverage<1) {
            Selected <- Posterior %>% arrange(CHR, BP.lead, desc(P.set)) %>%
                group_by(CHR, BP.lead) %>%
                mutate(P.set.sum = cumsum(P.set), Row.max=which.max(P.set.sum >= coverage), Row=row_number()) %>%
                filter(Row <= Row.max) %>%
                select(CHR, SNP.lead, BP.lead, SNP, BP, P.set, P.causal, P.set.sum, Causal) %>%
                ungroup()
        } else {
            Selected <- Posterior %>%
                select(CHR, SNP.lead, BP.lead, SNP, BP, P.set, P.causal, Causal)
        }
        # Summarise discoveries by clump
        Selected <- Selected %>% group_by(CHR, SNP.lead, BP.lead) %>%
            summarise(Causal=any(Causal), Size=n(), BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min) %>%
            ungroup() %>% arrange(CHR, BP.lead)
    } else {
        Selected <- tibble()
    }
    return(Selected)
}

Params.caviar <- expand.grid(Amplitude=amplitude.list, Fold=fold.list,
                             Clumping=clumping.list, Coverage=c(0.9)) %>%
    as_tibble() %>%
    mutate(Fold=as.character(Fold), Clumping=as.character(Clumping))

# Load posterior inclusion probabilities according to CAVIAR
Selections.caviar <- foreach(idx = 1:nrow(Params.caviar), .combine = 'rbind') %dopar% {
    amplitude <- Params.caviar$Amplitude[idx]
    fold <- Params.caviar$Fold[idx]
    clumping <- Params.caviar$Clumping[idx]
    coverage <- Params.caviar$Coverage[idx]
    cat(sprintf("Loading discoveries for fold %s, amplitude %d, clumping %s and refined with coverage %s...\n",
                fold, amplitude, clumping, coverage))
    Selections.row <- load.caviar(amplitude, fold, clumping, coverage=coverage) %>%
        mutate(Amplitude=amplitude, Fold=fold, Clumping=clumping, Method="LMM", Coverage=coverage)
    return(Selections.row)
}

# Summarise fine-mapping discoveries
Summary.caviar <- foreach(idx = 1:nrow(Params.caviar), .combine = 'rbind') %dopar% {
    fold <- Params.caviar$Fold[idx]
    clumping <- Params.caviar$Clumping[idx]
    amplitude <- Params.caviar$Amplitude[idx]
    coverage <- Params.caviar$Coverage[idx]
    Summary.row <- tibble(Clumping=clumping, Coverage=coverage, Amplitude=amplitude, Fold=fold)
    Selections.row <- Selections.caviar %>%
        inner_join(Summary.row, by = c("Clumping", "Amplitude", "Fold", "Coverage"))
    Summary.row$Discoveries <- nrow(Selections.row)
    Summary.row$FDP <- sum(!Selections.row$Causal)/max(1,Summary.row$Discoveries)
    Summary.row$Discoveries.true <- sum(Selections.row$Causal)
    Summary.row$Size <- median(Selections.row$Size,na.rm=T)
    Summary.row$BP.width <- as.integer(median(Selections.row$BP.width,na.rm=T))
    Summary.row$Signals <- sum(Variants.original$Causal)
    return(Summary.row)
}
Summary.caviar <- Summary.caviar %>%
    full_join(Architecture %>% select(Amplitude, H2), by = "Amplitude") %>%
    select(colnames(Summary.susie))

cat("Summary (fine-mapping with CAVIAR)\n")
Summary.caviar %>% filter(Clumping==clumping.list[1]) %>%
    group_by(Clumping, Amplitude, H2, Coverage) %>%
    summarise(Discoveries.true=mean(Discoveries.true), FDR=mean(FDP),
              Size=mean(Size), BP.width=mean(BP.width))

out.file <- sprintf("%s/simulations/summary/lmm_%s_caviar.txt", scratch, experiment.name)
Summary.caviar %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

################
## Make plots ##
################

if(FALSE) {

    Summary.susie %>%
        ggplot(aes(x=H2, y=Discoveries.true, color=Clumping)) + geom_point() + geom_line()

    Summary.caviar %>%
        ggplot(aes(x=H2, y=Discoveries.true, color=Clumping)) + geom_point() + geom_line()

    Summary.susie %>%
        ggplot(aes(x=H2, y=FDP, color=Clumping)) +
        geom_point() + geom_line() +
        ylim(0,1) +
        geom_hline(yintercept=0.1) +
        theme_bw()

        Summary.susie %>%
        ggplot(aes(x=H2, y=BP.width/1e6, color=Clumping)) +
        geom_point() + geom_line() +
        ylim(0,1) +
        theme_bw()

}
