library(tidyverse)
library(kableExtra)
source("../utils/utils_multires.R")

# Default arguments
experiment.name <- "eight"
experiment.numb <- "8"
fold <- "01"
conservative <- 0

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])
experiment.numb <- as.character(args[2])
fold <- as.character(args[3])
conservative <- as.integer(args[4])

## Global variables
scratch <- "/scratch/groups/candes/ukbiobank_tmp"
resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")
offset <- 1
seed <- 0

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

# Load list of causal variants
causal.file <- sprintf("%s/simulations_small/phenotypes/%s_causal.txt", scratch, experiment.name)
Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, BP)
Variants <- Variants %>%
    left_join(mutate(Causal,Causal=TRUE),by = c("CHR", "BP")) %>%
    mutate(Causal=replace_na(Causal,FALSE))

##########################
## Load test statistics ##
##########################

## Load knockoff statistics at each resolution
Stats <- lapply(resolution.list, function(resolution) {
  stats.file <- sprintf("%s/simulations_small/lasso/fold_%s/%s/stats_%s_%s.txt",
                        scratch, fold, resolution, experiment.name, experiment.numb)
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
    selected.idx <- which(data$groups[,m] %in% selected.groups)
    selected.variants <- filter(Variants, Resolution==resolution.list[m])[selected.idx,]
  }
})
Selections <- do.call("rbind", Selections)

## Summarise results
Results <- Selections %>%
  group_by(Resolution, CHR, Group) %>%
  summarise(BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n(), Causal=any(Causal)) %>%
  select(Resolution, CHR, BP.min, BP.max, BP.width, Size, Group, Causal) %>%
  arrange(Resolution, CHR, Group)

## Save results
discoveries.file <- sprintf("%s/simulations_small/lasso/fold_%s/multires/lasso_%s_%s_c%s.txt",
                            scratch, fold, experiment.name, experiment.numb, conservative)
Results %>% write_delim(discoveries.file, delim=" ")
cat(sprintf("List of multiresolution discoveries written to: %s\n", discoveries.file))
