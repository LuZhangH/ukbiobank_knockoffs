library(tidyverse)
library(kableExtra)
source("../utils/utils_multires.R")

## Input arguments
phenotype <- "bmi"
seed <- 0
conservative <- 0

## Input arguments
args <- commandArgs(trailingOnly=TRUE)
phenotype <- as.character(args[1])
seed <- as.integer(args[2])
conservative <- as.integer(args[3])

## Global variables
scratch <- "/scratch/groups/candes/ukbiobank_tmp"
resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")
offset <- 1

## Load list of variants
chr.list <- seq(1,22)
Params <- expand.grid(Resolution=resolution.list, CHR=chr.list) %>%
    as_tibble() %>% mutate(Resolution=as.character(Resolution), CHR=as.integer(CHR))
Variants <- lapply(1:nrow(Params), function(idx) {
    resolution <- Params$Resolution[idx]
    chr <- Params$CHR[idx]
    cat(sprintf("Loading list of variants at resolution %s for chromosome %d...\n", resolution, chr))
    ## Load list of variants
    grp.file <- sprintf("%s/clumping/%s/grp_chr%d.txt", scratch, resolution, chr)
    Variants.chr <- read_delim(grp.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group")
    Variants.chr$Resolution <- resolution
    ## Return results
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants)
Variants <- Variants %>% mutate(Resolution = factor(Resolution, levels=resolution.list, labels=resolution.list))

##########################
## Load test statistics ##
##########################

## Load knockoff statistics at each resolution
Stats <- lapply(resolution.list, function(resolution) {
  cat(sprintf("Loading discoveries for %s (seed %s) at resolution %s...\n", phenotype, seed, resolution))
  stats.file <- sprintf("%s/analysis/knockoffs/%s_%s_s%s_lasso_stats.txt",
                        scratch, phenotype, resolution, seed)
  Stats.res <- read_delim(stats.file, delim=" ", col_types=cols())
})
Stats <- do.call("rbind", Stats)
Stats <- Stats %>% mutate(Resolution = factor(Resolution, levels=resolution.list, labels=resolution.list))

## Remove stats with very large R2
Stats <- Stats %>% filter(is.na(R2)|R2<=0.99)

#############################
## Apply multilayer filter ##
#############################
data <- prepare_multires(Variants, Stats)

## Apply multilayer filter
if(conservative==0) {
  filter.output <- multilayer_knockoff_filter(data$W, data$groups, data$parents, 0.1, offset, "one-way")
} else {
  filter.output <- multilayer_knockoff_filter(data$W, data$groups, data$parents, 0.1/1.93, offset, "one-way")
}

## Extract list of discoveries from output
Selections <- lapply(1:length(data$W), function(m) {
  selected.groups <- filter.output$S_hats[[m]]
  if(length(selected.groups)>0) {
    Variants.res <- 
      selected.idx <- which(data$groups[,m] %in% selected.groups)
    selected.variants <- filter(Variants, Resolution==resolution.list[m])[selected.idx,]
  }
})
Selections <- do.call("rbind", Selections)

## Summarise results
Stats.tmp <- Stats %>% mutate(SNP=SNP.lead, BP=BP.lead) %>%
  select(Resolution, CHR, Group, SNP, BP, W)
Results <- Selections %>%
  group_by(Resolution, CHR, Group) %>%
  summarise(BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n()) %>%
  left_join(Stats.tmp, by = c("Resolution", "CHR", "Group")) %>%
  group_by(Resolution, CHR, Group) %>%
  summarise(BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n(),
            SNP.lead=SNP[which.max(W)], BP.lead=BP[which.max(W)], W=max(W)) %>%
  select(Resolution, CHR, SNP.lead, BP.lead, W, BP.min, BP.max, BP.width, Size, Group) %>%
  arrange(Resolution, desc(W))
    
## Save results
discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_multires_s%s_offset%s_c%s.txt",
                            scratch, phenotype, seed, offset, conservative)
Results %>% write_delim(discoveries.file, delim=" ")
cat(sprintf("List of multiresolution discoveries written to: %s\n", discoveries.file))

if(FALSE) {

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

  W.thres <- knockoff.threshold(filter(Stats, Resolution=="Radj2")$W, fdr=0.1, offset=1)

  Stats %>% filter(Resolution=="Radj2", W >= W.thres)

  Stats %>% filter(Resolution=="Radj2") %>% group_by(CHR,Group)
}
