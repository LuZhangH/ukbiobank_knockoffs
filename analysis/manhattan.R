.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))

source("utils_manhattan.R")

scratch <- "/scratch/PI/candes/ukbiobank_tmp"
resolution <- "Radj2"

# Load list of variants
Variants <- lapply(1:22, function(chr) {    # Load list of variants
    cat(sprintf("Loading list of variants on chromosome %d... ", chr))
    key.file <- sprintf("%s/knockoffs/%s_K50/ukb_gen_chr%d.key", scratch, resolution, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    # Load LD table
    ld.file <- sprintf("%s/knockoff_diagnostics/%s_K50/ukb_gen_chr%d.ld", scratch, resolution, chr)
    LD.chr <- read_table(ld.file, col_types=cols(), guess_max=21474836) %>%
        filter(BP_A==BP_B) %>%
        mutate(CHR=CHR_A, BP=BP_A) %>%
        select(CHR, BP, R2)
    # Combine list of variants with LD and MAF tables
    Variants.chr <- Variants.chr %>% left_join(LD.chr, by = c("CHR", "BP"))
    cat("done.\n")
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants)

# Specify phenotype
phenotype <- "respiratory"

# Load knockoff statistics
lasso.file <- sprintf("%s/analysis/knockoffs/%s_%s_lasso.txt", scratch, phenotype, resolution)
Lasso <- read_delim(lasso.file, delim=" ", col_types=cols())

# Compute the knockoff statistics
W.stats <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
}

Stats <- Lasso %>% select("Resolution", "CHR", "Group", "SNP", "BP", "Importance") %>%
    filter(Importance!=0) %>% 
    left_join(Variants, by = c("CHR", "Group", "SNP", "BP")) %>%
    group_by(Resolution, CHR, Group) %>%
    summarize(W = W.stats(Importance,Knockoff),
              Lead=which.max(Importance), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
              Size=n(), R2=max(R2)) %>%
    mutate(SNP.lead = gsub(".A", "", SNP.lead), SNP.lead = gsub(".B", "", SNP.lead)) %>%
    ungroup() %>%
    arrange(desc(abs(W))) %>%
    select(CHR, Group, SNP.lead, BP.lead, Size, W, R2, Resolution) %>%
    filter(W!=0)

# Save stats
stats.file <- sprintf("%s/analysis/knockoffs/%s_stats_%s.txt", scratch, phenotype, resolution)
Stats %>% write_delim(stats.file, delim=" ")
