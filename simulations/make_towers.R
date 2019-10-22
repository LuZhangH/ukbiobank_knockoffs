#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

experiment.name <- "seven"
amplitude <- 5

args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])
amplitude <- as.character(args[2])

library(tidyverse)
library(purrrlyr)
library(doParallel)
source("../utils/utils_clumping.R")
source("utils_outer.R")

K <- 50
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
#oak <- "/oak/stanford/groups/candes/ukbiobank_tmp"
chr.list <- seq(1, 22)
resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")

###########################
## Load list of variants ##
###########################

no_cores <- 1
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

##############################
## Load knockoff statistics ##
##############################

Params <- expand.grid(Amplitude=amplitude, Fold=1, Resolution=resolution.list) %>%
    as_tibble() %>%
    mutate(Resolution=as.character(Resolution))

Selections <- foreach(idx = 1:nrow(Params), .combine = 'rbind') %dopar% {
    fold <- Params$Fold[idx]
    amplitude <- Params$Amplitude[idx]
    resolution <- Params$Resolution[idx]
    cat(sprintf("Loading knockoff statistics for fold %s, amplitude %d at resolution %s...\n",
                fold, amplitude, resolution))
    # Load the lasso coefficients
    res.file <- sprintf("%s/simulations/lasso/fold_%s/%s/lasso_%s_%s.txt",
                        scratch, fold, resolution, experiment.name, amplitude)
    if(file.exists(res.file)){
        #Params.loaded <- rbind(Params.loaded, Params[idx,])
        Lasso.res <- read_delim(res.file, delim=" ", col_types=cols())
        Stats <- Lasso.res %>% select("CHR", "Group", "SNP", "BP", "Importance") %>%
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

        # Save stats file for the lowest resolution
        if(resolution=="Radj2") {
            stats.file <- sprintf("%s/simulations/summary/towers_%s_%s_stats.txt",
                                scratch, experiment.name, amplitude)
            Stats %>% write_delim(stats.file, delim=" ")
            cat(sprintf("Written stats for resolution %s on:\n %s\n", resolution, stats.file))            
        }
        
        # Apply the knockoff filter
        Selections.row <- Stats %>% filter(R2.max<0.99) %>% knockoff.filter(fdr=0.1, offset=1)
        Selections.row <- Selections.row %>%
            select(Resolution, CHR, SNP.lead, BP.lead, Group, W) %>%
            inner_join(Variants, by = c("Resolution", "CHR", "Group")) %>%
            group_by(Resolution, CHR, SNP.lead, BP.lead, Group) %>%
            summarise(W=mean(W), BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n()/2,
                      Causal=any(Causal)) %>%
            ungroup() %>%
            mutate(Fold=fold, Amplitude=amplitude, Method="Knockoffs")
    } else {
        warning(sprintf("File %s does not exist.\n", res.file))
        return(tibble())
    }
    return(Selections.row)
}

resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")
Selections <- Selections %>%
    mutate(Importance=W,
           Resolution=factor(Resolution, levels=resolution.list, labels=resolution.list))

##############################
## LMM p-values             ##
##############################

# Load LMM p-values
lmm.file <- sprintf("%s/simulations/bolt/fold_1/stats_%s_%s.txt", scratch, experiment.name, amplitude)
LMM <- read_delim(lmm.file, delim="\t", col_types=cols()) %>% as_tibble()
if("P_BOLT_LMM" %in% colnames(LMM)) {
    LMM <- LMM %>% mutate(P=P_BOLT_LMM)
} else {
    LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
}
LMM <- LMM %>% inner_join(Variants.original, by = c("SNP", "CHR", "BP"))

# Load clumped LMM results
clumping <- "0.00000005"
clumped.file <- sprintf("%s/simulations/bolt_clumped/fold_1/stats_%s_%s_%s.clumped",
                    scratch, experiment.name, amplitude, clumping)
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
        select(CHR, SNP.lead, BP.lead, P, SNP, BP)
    # Combine lists of lead and secondary SNPs
    Clumped <- rbind(Clumped.raw.lead, Clumped.raw.secondary) %>% as_tibble() %>%
        arrange(CHR, BP.lead, BP)
    # Cross-reference with complete list of variants
    Clumped <- Clumped %>%
        inner_join(select(Variants.original,CHR,SNP,BP,Causal), by = c("CHR", "SNP", "BP")) %>%
        select(-P)
}

##################
## Save results ##
##################

selections.file <- sprintf("%s/simulations/summary/towers_%s_%s_knockoffs.txt",
                    scratch, experiment.name, amplitude)
Selections %>% write_delim(selections.file, delim=" ")
cat(sprintf("Written knockoff results on:\n %s\n", selections.file))

lmm.file <- sprintf("%s/simulations/summary/towers_%s_%s_lmm.txt",
                    scratch, experiment.name, amplitude)
LMM %>% write_delim(lmm.file, delim=" ")
cat(sprintf("Written LMM results on:\n %s\n", lmm.file))

clumped.file <- sprintf("%s/simulations/summary/towers_%s_%s_clumped.txt",
                    scratch, experiment.name, amplitude)
Clumped %>% write_delim(clumped.file, delim=" ")
cat(sprintf("Written clumped LMM results on:\n %s\n", clumped.file))

###############
## Make plot ##
###############
# Pick a discovery and define a window around it
window.lead <- Selections %>% filter(Method=="Knockoffs", CHR==1, Resolution=="Radj100") %>%
    arrange(desc(Importance)) %>% head(1)
window.lead

window.width <- 5*1e5
window.chr <- window.lead$CHR
window.left <- window.lead$BP.lead - window.width
window.right <- window.lead$BP.lead + window.width
cat(sprintf("Selected window on chromosome %d between %.2f and %.2f Mb.\n",
            window.chr, window.left/1e6, window.right/1e6))

source("utils/tower_plots.R")
pp <- plot_sears_tower(window.chr, window.left, window.right, Selections, LMM, Clumped)
pp
