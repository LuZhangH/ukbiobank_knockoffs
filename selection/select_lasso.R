#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

# Input arguments
phenotype <- as.character(args[1])
clumping <- as.character(args[2])

library(tidyverse)
library(snpStats)
library(oem)
library(glmnet)
library(knockoff)
library(doParallel)
no_cores <- 10
cl <- makeCluster(no_cores, type="FORK", outfile="")
registerDoParallel(cl)

K <- 50
tmp.dir <- "/scratch/PI/candes/ukbiobank_tmp"
#chr.list <- seq(1,22)
chr.list <- seq(6,6)
n.PCs <- 5

# Load fam file
fam.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.fam", tmp.dir, clumping, K, 22)
Subjects <- read_delim(fam.file, delim=" ", col_types=cols(),
                       col_names = c("FID", "IID", "X1", "X2", "Sex", "X3"))

# Load list of variants
Variants <- tibble()
for(chr in chr.list) {
    cat(sprintf("Loading list of variants on chromosome %d... ", chr))
    key.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.key", tmp.dir, clumping, K, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants <- rbind(Variants, Variants.chr)
    cat("done.\n")
}

# Load p-values from marginal testing
Results <- tibble()
for(chr in chr.list) {
    cat(sprintf("Loading results from chromosome %d... ", chr))
    pvals.file <- sprintf("%s/association/%s/%s_K%d/ukb_chr%d.assoc",
                         tmp.dir, phenotype, clumping, K, chr)
    Results.chr <- read_table(pvals.file, col_types=cols(), guess_max=5000000)
    Results <- rbind(Results, Results.chr)
    cat("done.\n")
}
Results <- Results %>% select(CHR, SNP, BP, P) %>%
    inner_join(Variants, by=c("CHR", "SNP", "BP"))

# Symmetrize p-values for screening
Results.symmetric <- Results %>% group_by(CHR, BP, Group) %>% summarise(P=min(P)) %>%
    inner_join(select(Results, CHR, SNP, BP, Group, Knockoff), by=c("CHR", "BP", "Group")) %>%
    select(colnames(Results)) %>% ungroup()

# Filter groups based on marginal p-values
Groups.filtered <- Results.symmetric %>% group_by(CHR, Group) %>%
    summarise(BP = BP[which.min(replace_na(P,1))], P=min(P), Group.size=n()) %>%
    select(CHR, BP, Group, Group.size, P) %>%
    mutate(Importance = -log10(P)) %>%
    arrange(desc(Importance)) %>% head(300)

Variants.filtered <- Results.symmetric %>%
    inner_join(select(Groups.filtered, CHR, Group), by = c("CHR", "Group")) %>%
    group_by(CHR, Group, Knockoff) %>%
    arrange(P) %>% mutate(Rank=seq(length(P))) %>% ungroup() %>%
    arrange(Rank, P) %>% head(3000)

# Restore non-symmetrized p-values
Variants.filtered <- Variants.filtered %>% select(-c("Rank", "P")) %>%
    inner_join(select(Results, CHR, SNP, BP, P), by=c("CHR", "SNP", "BP"))

n.groups.filtered <- Variants.filtered %>% group_by(CHR, Group) %>% summarise(N=n()) %>% nrow
n.knockoffs.filtered <- Variants.filtered %>% filter(Knockoff==TRUE) %>% nrow
cat(sprintf("Number of variants that passed screening: %d\n", nrow(Variants.filtered)))
cat(sprintf("Number of knockoffs that passed screening: %d\n", n.knockoffs.filtered))
cat(sprintf("Number of groups that passed screening: %d\n", n.groups.filtered))

# Load phenotypes
pheno.file <- sprintf("%s/association/phenotypes/phenotypes.tab", tmp.dir)
Phenotypes <- read_delim(pheno.file, delim="\t", col_types=cols())
Phenotypes <- Subjects %>% select(FID, IID) %>% inner_join(Phenotypes, by = c("FID", "IID"))

# Extract response vector and matrix of covariates
y <- Phenotypes[[phenotype]]
X.cov <- Phenotypes %>% select(sex, paste("PC.",seq(1,n.PCs),sep="")) %>% as.matrix

# Impute missing values in the matrix of covariates
for(j in 1:ncol(X.cov)) {
    missing.rows <- which(is.na(X.cov[,j]))
    X.cov[missing.rows, j] <- mean(X.cov[,j], na.rm = TRUE)
}

# Load the genotypes for the variants that passed filtering
chr.list <- sort(unique(Variants.filtered$CHR))
X <- foreach(chr.idx=1:length(chr.list), .combine = 'cbind') %dopar% {
    chr <- chr.list[chr.idx]
    snplist.chr <- (Variants.filtered %>% filter(CHR==chr))$SNP
    cat(sprintf("Loading %d variants from chromosome %d... \n", length(snplist.chr), chr))
    basename.chr <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d", tmp.dir, clumping, K, chr)
    Dat <- read.plink(basename.chr, select.snps=snplist.chr)
    cat(sprintf("Loaded %d variants from chromosome %d. \n", length(snplist.chr), chr))
    X.chr <- as(Dat$genotypes, "numeric")    
}
gc()

# Fit the lasso
row.idx <- which(!is.na(y))

cat(sprintf("There are %d missing values in the response vector.\n", sum(is.na(y))))
cat(sprintf("Fitting the lasso:"))
cat(sprintf("   n = %d\n", length(row.idx)))
cat(sprintf("   p = %d SNPs + %d covariates\n", ncol(X), ncol(X.cov)))
variable.order <- c(seq(ncol(X.cov)),ncol(X.cov)+sample(ncol(X)))
penalty.factor <- c(rep(0,ncol(X.cov)),rep(1,ncol(X)))
X <- scale(cbind(X.cov,X))

lasso.fit <- xval.oem(X[row.idx,variable.order], y[row.idx], penalty="lasso",
                      penalty.factor=penalty.factor,
                      ncores=no_cores, nfolds=5,
                      nlambda=50, tol = 1e-05)
cat("done.\n")
summary(lasso.fit)
# Extract coefficients and cross-reference with list of variants
Results.lasso <- tibble(SNP=colnames(X)[variable.order],
                        Z=predict(lasso.fit, s="lambda.min", type="coefficients")[-1]) %>%
    left_join(Variants.filtered, by="SNP") %>%
    mutate(Importance=abs(Z)) %>% arrange(desc(Importance))

# Save results
s.knockoffs.file <- sprintf("%s/association/%s/%s_K%d/lasso.txt", tmp.dir, phenotype, clumping, K)
write_delim(Results.lasso, s.knockoffs.file, delim=" ")

# Compute the knockoff statistics
W.group.stats <- function(importance, knockoff) {
    Z.X  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    Z.Xk <- sum(importance[which(knockoff==TRUE)], na.rm=T)    
    w = Z.X-Z.Xk
}
Results.knockoffs <- Results.lasso %>% filter(! SNP %in% colnames(X.cov)) %>%
    group_by(CHR, Group) %>%
    summarize(W = W.group.stats(Importance,Knockoff), SNP=SNP[1], BP=BP[1], P=min(P)) %>%
    ungroup() %>%
    mutate(Sign.W=factor(sign(W),levels=c(-1,0,1))) %>%
    arrange(desc(abs(W))) %>%
    select(CHR, SNP, BP, Group, W, Sign.W, P)
Results.knockoffs %>% head(20)


# Apply the knockoff filter
W.thres <- knockoff.threshold(Results.knockoffs$W, offset=0)
Selected.knockoffs <- Results.knockoffs %>% filter(W >= W.thres)
cat(sprintf("Knockoffs: %d groups were selected.\n", nrow(Selected.knockoffs)))
