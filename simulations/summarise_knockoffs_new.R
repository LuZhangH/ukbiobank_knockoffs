#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

args <- commandArgs(trailingOnly=TRUE)

experiment.name <- "six"
experiment.name <- as.character(args[1])

library(tidyverse)
library(purrrlyr)
library(doParallel)
source("../utils/utils_clumping.R")

K <- 50
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1, 22)
resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")
#resolution.list <- c("Radj2", "Radj5", "Radj10")

amplitude.list <- seq(1,10)
fold.list <- seq(1)

###########################
## Load list of variants ##
###########################

no_cores <- 4
cl <- makeCluster(no_cores, type="FORK", outfile="")
registerDoParallel(cl)

Params <- expand.grid(Resolution=resolution.list, CHR=chr.list) %>%
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
    Variants.chr$Resolution <- resolution
    # Load LD table
    ld.file <- sprintf("%s/knockoff_diagnostics/%s_K%d/ukb_gen_chr%d.ld", scratch, resolution, K, chr)
    LD.chr <- read_table(ld.file, col_types=cols(), guess_max=Inf) %>%
        filter(BP_A==BP_B) %>%
        mutate(CHR=CHR_A, BP=BP_A) %>%
        select(CHR, BP, R2)
    # Load MAF table
    frq.file <- sprintf("%s/knockoff_diagnostics/%s_K50/ukb_gen_chr%d.frq", scratch, resolution, chr)
    MAF.chr <- read.table(frq.file, header=T) %>% as_tibble() %>% select(CHR, SNP, MAF)
    # Combine list of variants with LD and MAF tables
    Variants.chr <- Variants.chr %>% left_join(LD.chr, by = c("CHR", "BP")) %>%
        left_join(MAF.chr, by = c("CHR", "SNP"))
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
    filter(Resolution==resolution.list[1], Knockoff==FALSE) %>%
    mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP))

###############################
## Load genetic architecture ##
###############################
arch.file <- sprintf("%s/simulations/phenotypes/architecture.txt", scratch)
Architecture <- read_delim(arch.file, delim=" ", col_types=cols()) %>%
    filter(Name==experiment.name)

##########################################
## Load knockoff stats and apply filter ##
##########################################

knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset>1 | offset<0) {
    stop('Input offset must be between 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

knockoff.filter.p <- function(Stats, fdr=0.10, offset=1) {
  if(offset>1 | offset<0) {
    stop('Input offset must be between 0 or 1')
  }
  Stats <- Stats %>% arrange(P)
  ratio = sapply(seq(nrow(Stats)), function(j) {
      (offset + sum(Stats$W[1:j] < 0)) / max(1, sum(Stats$W[1:j] > 0))
  })
  is.ok <- which(ratio<=fdr)
  if(length(is.ok)==0) {
      Selected <- Stats %>% head(0)
  } else {
      Selected <- Stats %>% head(max(is.ok)) %>% filter(W>0)
  }
  return(Selected)
}

knockoff.filter <- function(Stats, fdr=0.1, offset=1) {
    W.thres <- knockoff.threshold(Stats$W, fdr=fdr, offset=offset)
    Selected <- Stats %>% filter(W >= W.thres)
    return(Selected)
}

# Compute the knockoff statistics
W.stats <- function(importance, knockoff) {
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
}

##########################
## Outer-node knockoffs ##
##########################

amplitude.list <- c(10)

Params <- expand.grid(Amplitude=amplitude.list, Fold=fold.list, Resolution=resolution.list) %>%
    as_tibble() %>%
    mutate(Resolution=as.character(Resolution))

Stats <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    fold <- Params$Fold[idx]
    amplitude <- Params$Amplitude[idx]
    resolution <- Params$Resolution[idx]
    cat(sprintf("Loading knockoff statistics for fold %s, amplitude %d at resolution %s...\n",
                fold, amplitude, resolution))

    # Load the lasso coefficients
    res.file <- sprintf("%s/simulations/lasso/fold_%s/%s/lasso_%s_%s.txt",
                        scratch, fold, resolution, experiment.name, amplitude)
    Lasso.res <- read_delim(res.file, delim=" ", col_types=cols())

    Stats.row <- Lasso.res %>% select("CHR", "Group", "SNP", "BP", "Importance") %>%
        left_join(Variants, by = c("CHR", "Group", "SNP", "BP")) %>%
        mutate(R2=replace_na(R2,0), MAF=replace_na(MAF,0)) %>%
        group_by(Resolution, CHR, Group) %>%
        summarize(W = W.stats(Importance,Knockoff),
                  SNP.lead=SNP[1], BP.lead=BP[1],
                  MAF=min(MAF),
                  Causal=any(Causal),
                  R2.max=max(R2,na.rm=T), R2=mean(R2,na.rm=T),
                  Size=n()/2) %>%
        ungroup() %>%
        arrange(desc(abs(W))) %>%
        select(CHR, SNP.lead, BP.lead, Group, Size, Causal, W, MAF, R2.max, Resolution) %>%
        filter(W!=0)

    return(Stats.row)
}

Selections <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    resolution <- Params$Resolution[idx]
    cat(sprintf("Applyting knockoff filter for resolution %s...\n", resolution))

    Stats.res <- Stats %>% filter(Resolution==resolution)

    # Apply the knockoff filter
    Selections.row <- Stats.res %>% filter(R2.max<0.99) %>% knockoff.filter(fdr=0.1, offset=1)

    Selections.row <- Selections.row %>%
        select(Resolution, CHR, SNP.lead, BP.lead, Group, W) %>%
        inner_join(Variants, by = c("Resolution", "CHR", "Group")) %>%
        group_by(Resolution, CHR, SNP.lead, BP.lead, Group) %>%
        summarise(W=mean(W), BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n()/2,
                  Causal=any(Causal)) %>%
        ungroup()

    return(Selections.row)
}

source("utils_outer.R")

# Apply outer node filter
Params.outer <- expand.grid(Fold=fold.list, Amplitude=amplitude.list) %>%
    as_tibble() %>% mutate(Fold=as.character(Fold))

Selections.outer <- foreach(idx = 1:nrow(Params.outer), .combine = 'rbind') %dopar% {
    fold <- Params$Fold[idx]
    amplitude <- Params$Amplitude[idx]
    cat(sprintf("Resolution filter for amplitude %d and fold %s ... \n", amplitude, fold))
    Selections.outer.row <- Selections %>% mutate(Resolution=parse_number(Resolution)) %>%
        filter(Amplitude==amplitude, Fold==fold) %>% outer_filter()
}

Selections.outer %>% summarise(N=n(), True=sum(Causal), False=sum(!Causal), FDP=mean(!Causal))

# Compute FDP after outer-node filter
Selections.outer %>% group_by(Resolution) %>%
    summarise(N=n(), True=sum(Causal), False=sum(!Causal), FDP=mean(!Causal))

# Estimate FDP after outer filter
Errors <- estimate_fdp_outer(Stats)

Selections.outer %>% group_by(Resolution) %>%
    summarise(N=n(), True=sum(Causal), False=sum(!Causal), FDP=mean(!Causal)) %>%
    inner_join(Errors %>% mutate(Resolution=parse_number(Resolution))) %>%
    mutate(FDP.hat=False.hat/N)

Selections %>% group_by(Resolution) %>%
    summarise(N=n(), True=sum(Causal), False=sum(!Causal), FDP=mean(!Causal))

stopCluster(cl)
