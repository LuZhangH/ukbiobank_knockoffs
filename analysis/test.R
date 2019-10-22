#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

# Install packages bigsnpr and bigstatsr
# devtools::install_github("privefl/bigstatsr")
# devtools::install_github("privefl/bigsnpr")

# Documentation here:
# https://privefl.github.io/bigsnpr/reference/index.html
# https://privefl.github.io/bigstatsr/reference/big_spLinReg.html

# Load packages
library(tidyverse)

set.seed(1)
scratch <- "/scratch/PI/candes/ukbiobank_tmp"

chr.list <- seq(1,22)

###########################
## Load list of variants ##
###########################

# Load list of variants
Variants <- tibble()
for(chr in chr.list) {
    cat(sprintf("Loading list of variants on chromosome %d... ", chr))
    key.file <- sprintf("%s/knockoffs/%s_K50/ukb_gen_chr%d.key", scratch, resolution, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants <- rbind(Variants, Variants.chr)
    cat("done.\n")
}

# Compute the knockoff statistics
W.stats <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
}

###############################
## Apply the knockoff filter ##
###############################
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
knockoff.filter <- function(Stats, fdr=0.1, offset=1) {
    W.thres <- knockoff.threshold(Stats$W, fdr=fdr, offset=offset)
    Selected <- Stats %>% filter(W >= W.thres)
    return(Selected)
}

resolution <- "Radj2"
phenotype <- "cvd"
phenotype.list <- c("height", "bmi","sbp","platelet","cvd","diabetes","hypothyroidism","respiratory","glaucoma")
phenotype.list <- c("height")

for( phenotype in phenotype.list) {
    
    # Load knockoff stats
    stats.file <- sprintf("%s/analysis/knockoffs/%s_%s_lasso.txt", scratch, phenotype, resolution)
    Lasso.res <- read_delim(stats.file, delim=" ") %>% mutate(Z=Importance) %>% select(-Importance)

    # 
    Stats <- Lasso.res %>% select("Resolution", "CHR", "Group", "SNP", "BP", "Z") %>%
        filter(Z!=0) %>% 
        left_join(Variants, by = c("CHR", "Group", "SNP", "BP")) %>%
        group_by(Resolution, CHR, Group) %>%
        summarize(W = W.stats(abs(Z),Knockoff),
                  Lead=which.max(abs(Z)),
                  SNP.lead=SNP[Lead], BP.lead=BP[Lead],
                  Size=n()) %>%
        mutate(SNP.lead = gsub(".A", "", SNP.lead), SNP.lead = gsub(".B", "", SNP.lead)) %>%
        ungroup() %>%
        arrange(desc(abs(W))) %>%
        select(CHR, Group, SNP.lead, BP.lead, Size, W, Resolution) %>%
        filter(W!=0)
    cat("done.\n")
    # Show a preview of the knockoff stats
    cat("Knockoff statistics:\n")
    Stats

    Selections <- Stats %>% knockoff.filter(fdr=0.1, offset=1)

    Selections <- Selections %>%
        select(Resolution, CHR, Group, SNP.lead, BP.lead, W) %>% inner_join(Variants, by = c("CHR", "Group")) %>%
        group_by(Resolution, CHR, Group, SNP.lead, BP.lead) %>%
        summarise(W=mean(W), BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n()/2,) %>%
        ungroup() %>%
        mutate(Method="Knockoffs") %>%
        arrange(desc(W))

    # Give a preview of the selections
    cat("Knockoff selections:\n")
    Selections %>% print()

    # Save list of discoveries
    discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_%s.txt", scratch, phenotype, resolution)
    Selections %>% select(CHR, SNP.lead, BP.lead, W, BP.min, BP.max, BP.width, Size, Group) %>%
        write_delim(discoveries.file, delim=" ")
    cat(sprintf("Saved list of %d discoveries to %s\n", nrow(Selections), discoveries.file))

}
