library(tidyverse)
library(gridExtra)
library(knockoff)
library(qqman)
library(scales)

phenotype <- "rd"
K <- 50
tmp.dir <- "/scratch/PI/candes/ukbiobank_tmp"
out.dir <- "/scratch/PI/candes/ukbiobank_tmp/results"
offset <- 1

load.results <- function(resolution) {
    # Load original list of variants
    Variants <- tibble()
    for(chr in seq(1,22)) {
        key.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.key", tmp.dir, resolution, K, chr)
        Variants.chr <- read_delim(key.file, delim=" ", col_types=cols(), progress=FALSE) %>%
            mutate(CHR=Chr, SNP=Variant, BP=Position) %>%
            select(CHR,SNP,BP,Group,Knockoff)
        corr.file <- sprintf("%s/augmented_data/%s_K%d/ukb_gen_chr%d.ld", tmp.dir, resolution, K, chr)
        Corr.chr <- read_table(corr.file, col_types=cols(), guess_max=1000000)
        Corr.chr <- Corr.chr %>% filter(BP_A==BP_B) %>% mutate(CHR=CHR_A, SNP=SNP_A, BP=BP_A) %>%
            select(CHR, BP, R)
        Variants.chr <- Variants.chr %>% left_join(Corr.chr, by=c("CHR", "BP"))
        Variants <- rbind(Variants, Variants.chr)
    }

    # Load lasso importance measures
    res.file <- sprintf("%s/association/%s/%s_K%d/lasso.txt", tmp.dir, phenotype, resolution, K)
    Lasso <- read_delim(res.file, delim=" ", col_types=cols(), guess_max=1000000)

    # Cross-reference with list of variants
    Lasso <- Lasso %>% inner_join(Variants, by = c("SNP", "CHR", "BP", "Group", "Knockoff")) %>%
        select(colnames(Lasso), everything())

    # Compute knockoff statistics
    W.stats <- function(importance, knockoff) {
        z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
        zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
        w = z-zk
    }
    Results.knockoffs <- Lasso %>% filter(!is.na(CHR)) %>%
        group_by(CHR, Group) %>%
        summarize(W = W.stats(Importance,Knockoff), P = min(P),
                  Sign.W = factor(sign(W), levels=c(-1,0,1)), R=mean(R)) %>%
        ungroup() %>%
        arrange(desc(abs(W))) %>%
        select(CHR, Group, W, Sign.W, P, R) %>%
        mutate(Resolution = resolution)

    #Results.knockoffs <- Results.knockoffs %>%
    #    mutate(W=W/(1+(abs(R)>0.95))) %>% arrange(desc(abs(W)))

    # Apply the knockoff filter
    W.thres <- knockoff.threshold(Results.knockoffs$W, offset=0)
    Selected.knockoffs <- Results.knockoffs %>% filter(W >= W.thres)

    # Cross reference discoveries with list of variants
    Discoveries <- Variants %>% filter(Knockoff==FALSE) %>%
        inner_join(Selected.knockoffs, by=c("CHR", "Group")) %>%
        group_by(CHR, Group, Resolution) %>%
        summarize(Size=length(BP),
                  BP.min=as.integer(min(BP)), BP.max=as.integer(max(BP)), BP=as.integer(round(mean(BP))),
                  SNP = SNP[which.min(P)], P=min(P), W=mean(W)) %>%
        arrange(desc(W)) %>%
        ungroup() %>%
        select(CHR, SNP, BP, BP.min, BP.max, Group, Size, P, W, Resolution)

    # Count groups
    Groups <- Variants %>% filter(Knockoff==FALSE) %>%
        group_by(CHR, Group) %>%
        summarise(Size=n()) %>%
        mutate(Resolution = resolution) %>%
        ungroup()

    # Combine and return output
    output <- c()
    output$Discoveries <- Discoveries
    output$Groups <- Groups
    return(output)
}

resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50")
Discoveries <- tibble()
Groups <- tibble()
for(resolution in resolution.list) {
    cat(sprintf("Loading results at resolution %s... ", resolution))
    results <- load.results(resolution)
    Discoveries <- rbind(Discoveries, results$Discoveries)
    Groups <- rbind(Groups, results$Groups)
    cat("done.\n")
}

# Save results
out.file <- sprintf("%s/%s.txt", out.dir, phenotype)
Discoveries %>% arrange(CHR, BP, Resolution) %>%
    write_delim(out.file, delim=" ")
out.file <- sprintf("%s/%s_groups.txt", out.dir, phenotype)
Groups %>% arrange(Resolution, CHR, Group) %>%
    write_delim(out.file, delim=" ")

# Count discoveries at different resolutions
cat("Discoveries:\n")
Discoveries %>%
    mutate(Resolution = as.numeric(str_extract(Resolution,"[[:digit:]]+"))) %>%
    group_by(Resolution) %>% summarise(N = n())
cat("Groups:\n")
Groups %>% group_by(Resolution) %>%
    summarise(N=n(), Size=mean(Size))

# Make zoomed-in plot
chr <- 6
bp.center <- filter(Discoveries, CHR==chr, Resolution=="Radj50")$BP[1]
bp.width  <- 20e6
Discoveries %>%
    mutate(Importance = W / sum(Discoveries$W[Discoveries$Group==Group])) %>%
    filter(CHR==chr) %>%
    filter(abs(BP-bp.center) < bp.width) %>%
    mutate(Resolution = factor(Resolution, levels=resolution.list, ordered=TRUE)) %>%
    mutate(Height = as.numeric(Resolution)) %>%
    ggplot(aes(color=Resolution)) +
    geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.45, ymax=Height+0.45, alpha=Importance), fill="gray") +
    scale_y_discrete(labels=c("1" = "Radj2", "2" = "Radj5", "3" = "Radj10", "4"="Radj25", "5" = "Radj50"),
                     limits=c("1","2","3","4","5")) +
    scale_fill_manual(values=resolution.list) +
    theme_bw() + theme(legend.position="bottom")

# Manhattan plot
Discoveries %>%
    filter(Resolution == "Radj50") %>%
    manhattan(chr="CHR", bp="BP", snp="SNP", p="W", logp=FALSE, genomewideline="FALSE", ylab="Importance")
