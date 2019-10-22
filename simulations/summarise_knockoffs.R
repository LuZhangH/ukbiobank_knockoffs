#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

args <- commandArgs(trailingOnly=TRUE)

experiment.name <- "six"
experiment.name <- as.character(args[1])

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

###########################
## Load list of variants ##
###########################

no_cores <- 8
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

Params <- expand.grid(Amplitude=amplitude.list, Fold=fold.list, Resolution=resolution.list) %>%
    as_tibble() %>%
    mutate(Resolution=as.character(Resolution))
#Params.loaded <- Params %>% head(0)

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
    summarise(Discoveries=sum(!is.na(SNP.lead)), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              BP.width=mean(BP.width), Size=mean(Size)) %>%
    mutate(Method="Knockoffs")

cat("Summary of distinct discoveries with Knockoffs:\n")
Distinct.summary %>% group_by(Amplitude, H2, Resolution) %>%
    summarise(Discoveries.true=mean(Discoveries.true), FDP=mean(FDP)) %>%
    filter(Resolution=="Radj2") %>% print()

out.file <- sprintf("%s/simulations/summary/knockoffs_%s_distinct.txt", scratch, experiment.name)
Distinct.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

Selections.loci.consolidated <- Selections.loci %>%
    mutate(Importance=W) %>%
    group_by(Amplitude, H2, Fold, Resolution) %>%
    consolidate_clumps(gap=1e5)
    
Distinct.consolidated.summary <- Selections.loci.consolidated %>%
    group_by(Amplitude, H2, Fold, Resolution) %>%
    summarise(Discoveries=sum(!is.na(SNP.lead)), Discoveries.true=sum(Causal,na.rm=T), FDP=mean(!Causal),
              BP.width=mean(BP.width), Size=mean(Size)) %>%
    mutate(Method="Knockoffs")

cat("Summary of distinct consolidated discoveries with knockoffs:\n")
Distinct.consolidated.summary %>% print()

out.file <- sprintf("%s/simulations/summary/knockoffs_%s_distinct_consolidated.txt", scratch, experiment.name)
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

out.file <- sprintf("%s/simulations/summary/knockoffs_%s_loci.txt", scratch, experiment.name)
Cluster.summary %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to %s \n", out.file))

stopCluster(cl)
