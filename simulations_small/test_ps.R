#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

library(tidyverse)
library(bigsnpr)

scratch <- "/scratch/PI/candes/ukbiobank_tmp/"

experiment <- "sevenps"
resolution <- "Radj2"
fold <- "01"
chr <- 1

###############################
## Load genetic architecture ##
###############################
arch.file <- sprintf("%s/simulations_small/phenotypes/architecture.txt", scratch)
Architecture <- read_delim(arch.file, delim=" ", col_types=cols()) %>%
    filter(Name==experiment)

###########################
## Load list of variants ##
###########################

chr.list <- seq(1,22)
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

# Load list of causal variants
causal.file <- sprintf("%s/simulations_small/phenotypes/%s_causal.txt", scratch, experiment)
Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, BP)
Variants <- Variants %>%
    left_join(mutate(Causal,Causal=TRUE),by = c("CHR", "BP")) %>%
    mutate(Causal=replace_na(Causal,FALSE)) %>%
    mutate(Causal=ifelse(Knockoff, FALSE, Causal))

####################
## Load genotypes ##
####################

tmp.file <- sprintf("%s/simulations_small/lasso_data/ukb_gen_%s_fold_%s", scratch, resolution, fold)
rds.file <- sprintf("%s.rds", tmp.file)
if(file.exists(rds.file)){
    cat(sprintf("Found FBM in %s.\n", rds.file))
} else {
    cat(sprintf("Could not find FBM in %s. Starting conversion. \n", rds.file))
    bed.file <- sprintf("%s/simulations_small/lasso_data/ukb_gen_%s_fold_%s.bed", scratch, resolution, fold)
    bk.file <- sprintf("%s.bk", tmp.file)
    if(file.exists(bk.file)) {
        file.remove(bk.file)
    }
    cat(sprintf("Temporary file for FBM: %s\n", bk.file))
    cat("Reading genotype file... ")
    snp_readBed(bed.file, backingfile = tmp.file)
    cat("done.\n")
}

# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach(rds.file)

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
SNP <- obj.bigSNP$map$marker.ID

##############
## Load PCs ##
##############

# Load fam file
fam.file <- sprintf("%s/knockoffs/Radj10_K50/ukb_gen_chr%d.fam", scratch, 22)
Subjects <- read_delim(fam.file, delim=" ", col_types=cols(),
                       col_names = c("FID", "IID", "X1", "X2", "Sex", "X3"))

# Load response values
pheno.file <- sprintf("%s/phenotypes/phenotypes_qc.tab", scratch)
Confounders <- read_delim(pheno.file, delim="\t", col_types=cols())
Confounders <- Confounders %>% right_join(Subjects, by=c("FID", "IID"))
# Make sure that the rows of the genotypes match the rows of the phenotypes
Confounders <- Confounders %>%
    right_join(transmute(obj.bigSNP$fam, FID=family.ID, IID=sample.ID), by = c("FID", "IID"))
# Extract confounder matrix
Z <- Confounders %>% select(starts_with("PC.")) %>% as.matrix() %>% scale()
# Impute missing values
Z[is.na(Z)] <- 0

# Cluster the PCs
k.list <- seq(10)
wss <- sapply(k.list, function(k) {
    cat(sprintf("Clustering with %d groups...\n", k))
    kmeans(Z, k, nstart=5)$tot.withinss
})

plot(k.list, wss, type="b",
     pch = 19, frame = FALSE, xlab="Number of clusters K", ylab="Total within-clusters sum of squares")

pcCluster <- kmeans(Z, 6)
table(pcCluster$cluster)

###################
## Load p-values ##
###################

qqplot.pvals <- function(pvals, pvals.lmm) {
    observed <- sort(pvals)
    lobs <- -(log10(observed))
    observed.lmm <- sort(pvals.lmm)
    lobs.lmm <- -(log10(observed.lmm))
    expected <- c(1:length(observed))
    lexp <- -(log10(expected / (length(expected)+1)))
    dat <- tibble(Marginal=lobs, LMM=lobs.lmm, expected=lexp) %>%
        gather(Marginal, LMM, key = "method", value = "nlogp")
    pp <- dat %>%
        ggplot(aes(x=expected, y=nlogp, color=method)) +
        geom_point(alpha=0.2) +
        geom_abline(intercept = 0, slope = 1) +
        xlim(0,10) + ylim(0,10) +
        xlab("Expected") + ylab("Observed") +
        theme_bw() +
        theme(legend.position = c(0.8, 0.2))
    return(pp)
}

# Load marginal pvalues
amplitude <- 5
pvals.file <- sprintf("%s/simulations_small/pvalues/fold_%s/%s_%s_chr%s.assoc.linear",
                      scratch, fold, experiment, amplitude, chr)
Pvals <- read_table(pvals.file, col_types=cols())

# Load LMM pvalues
lmm.file <- sprintf("%s/simulations_small/bolt/fold_%s/stats_%s_%s.txt",
                     scratch, fold, experiment, amplitude)
LMM.raw <- read_delim(lmm.file, delim="\t", col_types=cols())
LMM <- LMM.raw %>% mutate(P=P_BOLT_LMM) %>% filter(CHR==chr)

p.qqplot <- qqplot.pvals(Pvals$P, LMM$P)
plot.file <- sprintf("figures/qqplot_a%s.pdf", amplitude)
ggsave(plot.file, plot=p.qqplot, width=4, height=4)

#######################
## Summarize results ##
#######################

amplitude <- 1

# Load knockoff locus results
resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj100")
knock.res.file <- sprintf("%s/simulations_small/summary/knockoffs_sevenps_loci.txt", scratch)
Results.knock <- read_delim(knock.res.file, delim=" ", col_types=cols()) %>%
    mutate(Resolution = factor(Resolution, levels=resolution.list, labels=resolution.list))

# Load lmm locus results
lmm.res.file <- sprintf("%s/simulations_small/summary/lmm_sevenps_loci.txt", scratch)
Results.lmm <- read_delim(lmm.res.file, delim=" ", col_types=cols())

# Combine and plot results
Results <- rbind(
    Results.knock %>% filter(Resolution %in% c("Radj2", "Radj5", "Radj20", "Radj100")),
    Results.lmm %>% filter(Clumping %in% c("0.00000005", "0.00005")) %>% mutate(Resolution=paste("LMM (", Clumping, ")",sep="")) %>% select(-Clumping))

p.fdp <- Results %>%
    ggplot(aes(x=H2, y=FDP, color=Resolution, linetype=Method)) +
    geom_point() + geom_line() +
    geom_abline(intercept=0.1, slope=0) +
    ylim(0,1) +
    theme_bw() +
    theme(legend.position = "bottom") + 
    guides(color=guide_legend(ncol=1), linetype=guide_legend(ncol=1))
plot.file <- sprintf("figures/fdp.pdf", amplitude)
ggsave(plot.file, plot=p.fdp, width=5, height=6)

p.pow <- Results %>%
    ggplot(aes(x=H2, y=Power, color=Resolution, linetype=Method)) +
    geom_point() + geom_line() +
    ylim(0,1) +
    theme_bw() +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(ncol=1), linetype=guide_legend(ncol=1))
plot.file <- sprintf("figures/pow.pdf", amplitude)
ggsave(plot.file, plot=p.pow, width=5, height=6)

#######################################################################
## Proportion of variance in the PCs explained by selected variables ##
#######################################################################

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


# Load variables selected with knockoffs
amplitude <- 5
stats.file <- sprintf("%s/simulations_small/lasso/fold_%s/%s/lasso_%s_%s.txt",
                      scratch, fold, resolution, experiment, amplitude)
Lasso <- read_delim(stats.file, delim=" ", col_types=cols())

Stats <- Lasso %>% select("CHR", "Group", "SNP", "BP", "Importance") %>%
    left_join(Variants, by = c("CHR", "Group", "SNP", "BP")) %>%
    group_by(CHR, Group) %>%
    summarize(W = W.stats(Importance,Knockoff),
              SNP.lead=SNP[1], BP.lead=BP[1],
              Causal=any(Causal),
              Size=n()/2) %>%
    ungroup() %>%
    arrange(desc(abs(W))) %>%
    select(CHR, SNP.lead, BP.lead, Group, Size, Causal, W) %>%
    filter(W!=0) %>%
    inner_join(Variants %>% mutate(SNP.lead=SNP) %>% select(SNP.lead, Knockoff))

Selections <- Stats %>% knockoff.filter(fdr=0.1, offset=1)

# Compute PC variance explained by falsely selected variants
Selections.false <- Selections %>% filter(Causal==FALSE)
xf.idx <- which(SNP %in% Selections.false$SNP.lead)
Xf <- G[,xf.idx]
data <- cbind(Z, Xf) %>% as_tibble()
colnames(data) <- c(paste("PC", seq(ncol(Z)), sep="."), Selections.false$SNP.lead)
lm.pc <- lm(PC.1~.-PC.2-PC.3-PC.4-PC.5, data)
summary(lm.pc)

# Compute PC variance explained by selected variants
xf.idx <- which(SNP %in% Selections$SNP.lead)
Xf <- G[,xf.idx]
data <- cbind(Z, Xf) %>% as_tibble()
colnames(data) <- c(paste("PC", seq(ncol(Z)), sep="."), Selections.false$SNP.lead)
lm.pc <- lm(PC.1~.-PC.2-PC.3-PC.4-PC.5, data)
summary(lm.pc)
