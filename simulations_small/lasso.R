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
library(devtools)
library(bigsnpr)

# Default arguments
experiment.name <- "zero"
experiment.numb <- "10"
fold <- "01"
resolution <- "Radj100"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])
experiment.numb <- as.character(args[2])
fold <- as.character(args[3])
resolution <- as.character(args[4])

# Other parameters
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1,22)
ncores <- 5
fit.lasso <- FALSE
seed <- 0

###########################
## Load list of variants ##
###########################

# Load list of variants
Variants <- tibble()
for(chr in chr.list) {
    cat(sprintf("Loading list of variants on chromosome %d... ", chr))
    key.file <- sprintf("%s/knockoffs/%s_K50_s%s/ukb_gen_chr%d.key", scratch, resolution, seed, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    # Load LD table
    ld.file <- sprintf("%s/knockoff_diagnostics/%s_K50_s%s/ukb_gen_chr%d.ld", scratch, resolution, seed, chr)
    LD.chr <- read_table(ld.file, col_types=cols(), guess_max=Inf) %>%
        filter(BP_A==BP_B) %>%
        mutate(CHR=CHR_A, BP=BP_A) %>%
        select(CHR, BP, R2)
    # Combine list of variants with LD and MAF tables
    Variants.chr <- Variants.chr %>% left_join(LD.chr, by = c("CHR", "BP"))
    Variants <- rbind(Variants, Variants.chr)
    cat("done.\n")
}

# Load list of causal variants
causal.file <- sprintf("%s/simulations_small/phenotypes/%s_causal.txt", scratch, experiment.name)
Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, BP)
Variants <- Variants %>%
    left_join(mutate(Causal,Causal=TRUE),by = c("CHR", "BP")) %>%
    mutate(Causal=replace_na(Causal,FALSE)) %>%
    mutate(Causal=ifelse(Knockoff, FALSE, Causal))

####################
## Load genotypes ##
####################

if(fit.lasso) {

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

  ## Attach the "bigSNP" object in R session
  obj.bigSNP <- snp_attach(rds.file)

  ## Get aliases for useful slots
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  SNP <- obj.bigSNP$map$marker.ID
  
}

# Compute scaling factor for the genotypes
scale.file <- sprintf("%s/simulations_small/lasso_data/ukb_gen_%s_fold_%s_scale.txt", scratch, resolution, fold)
if(!file.exists(scale.file)){
    ptm <- proc.time() # Start the clock!
    cat("Computing scaling factors for all variants... ")
    scaler <- big_scale()
    G.scale <- scaler(G)
    scaling.factors <- G.scale$scale
    cat("done.\n")
    proc.time() - ptm # Stop the clock
    # Save the scaling factors to file
    write.table(scaling.factors, scale.file, sep = "\t", row.names=F, col.names=F)
    cat(sprintf("Saved scaling factors to %s\n", scale.file))
} else {
    cat(sprintf("Loading scaling factors from %s... ", scale.file))
    scaling.factors <- read.table(scale.file)$V1
    cat("done.\n")
}

#####################
## Load phenotypes ##
#####################

if(fit.lasso) {

  cat("Reading phenotype file... ")

  ## Load list of subjects
  Subjects <- obj.bigSNP$fam %>% as_tibble() %>% transmute(FID=family.ID, IID=sample.ID)

  ## Load response values
  pheno.file <- sprintf("%s/simulations_small/phenotypes/%s_phenotypes_%s.tab", scratch, experiment.name, fold)
  Phenotypes <- read_delim(pheno.file, delim="\t", col_types=cols())
  keep.rows <- which(Phenotypes$FID %in% Subjects$FID)
  Phenotypes <- Phenotypes %>% right_join(Subjects, by=c("FID", "IID"))

  cat("done.\n")
  
}


###################
## Fit the lasso ##
###################

if(fit.lasso) {

  dfmax <- 20000

  ## Extract phenotype
  phenotype <- sprintf("Y_a%s_1", experiment.numb)
  y <- Phenotypes[[phenotype]]

  ## Check whether the phenotype is continuous or binary
  if(all(y %in% c(1,2))) {
    y <- factor(y, levels=c(1,2), labels=c(0,1))
    y <- as.numeric(levels(y))[y]
    y.binary <- TRUE
    cat(sprintf("The phenotype is binary: %d cases and %d controls.\n", sum(y==1), sum(y==0)))
  } else {
    y.binary <- FALSE
  }
  
  ## Fit the lasso
  if(y.binary) {
    ptm <- proc.time() # Start the clock!
    cat(sprintf("Fitting sparse logistic regression with %d observations and %d variables... ",
                length(y), ncol(G)))
    lasso.fit <- big_spLogReg(G, y01.train=y, dfmax=dfmax, ncores=ncores)
    cat("done.\n")
    proc.time() - ptm # Stop the clock
  } else {
    ptm <- proc.time() # Start the clock!
    cat(sprintf("Fitting sparse linear regression with %d observations and %d variables... ",
                length(y), ncol(G)))
    lasso.fit <- big_spLinReg(G, y.train=y, dfmax=dfmax, ncores=ncores)
    cat("done.\n")
    proc.time() - ptm # Stop the clock
  }

  ## Extract beta from each fold and combine them
  cat("Extracting importance measures... ")
  beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)
  ## Undo scaling of lasso coefficients
  beta <- beta * scaling.factors
  Beta <- cbind(tibble(CHR=Variants$CHR,
                       SNP=Variants$SNP, BP=Variants$BP),
                as_tibble(beta)) %>% as_tibble()
  colnames(Beta) <- c("CHR", "SNP", "BP", paste("K", seq(ncol(beta)),sep=""))
  Beta <- Beta %>%
    mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
           Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) #%>%
  ##    filter(Nonzero==10)
  ## Extract the estimated coefficients
  Lasso.res <- Beta %>% mutate(Resolution=resolution) %>%
    inner_join(Variants, by = c("CHR", "SNP", "BP")) %>%
    filter(Z!=0) %>%
    arrange(desc(Z))
  cat("done.\n")

  ## Print summary of fit
  Lasso.res %>% summarise(Nonzero=sum(Z!=0), Causal=sum(Causal))

  ## Save lasso estimates
  out.file <- sprintf("%s/simulations_small/lasso/fold_%s/%s/lasso_%s_%s.txt",
                      scratch, fold, resolution, experiment.name, experiment.numb)
  Lasso.res %>% mutate(Importance=abs(Z)) %>%
    select(Resolution, CHR, SNP, BP, Group, Importance) %>%
    write_delim(out.file, delim=" ")
  cat(sprintf("Results saved to %s.\n", out.file))

} else {
  out.file <- sprintf("%s/simulations_small/lasso/fold_%s/%s/lasso_%s_%s.txt",
                      scratch, fold, resolution, experiment.name, experiment.numb)
  Lasso.res <- read_delim(out.file, delim=" ", col_types=cols()) %>%
    mutate(Z=Importance)
}

##################################
## Assemble knockoff statistics ##
##################################
cat("Computing knockoff statistics... ")

# Compute the knockoff statistics
W.stats <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
}

Stats <- Lasso.res %>% select("Resolution", "CHR", "Group", "SNP", "BP", "Z") %>%
    filter(Z!=0) %>%
    left_join(Variants, by = c("CHR", "Group", "SNP", "BP")) %>%
    group_by(Resolution, CHR, Group) %>%
    summarize(W = W.stats(abs(Z),Knockoff),
              SNP.lead=SNP[1], BP.lead=BP[1],
              Causal=any(Causal),
              R2=max(R2, na.rm=T),
              Size=n()) %>%
    ungroup() %>%
    arrange(desc(abs(W))) %>%
    select(CHR, Group, Size, Causal, W, R2, Resolution) %>%
    filter(W!=0)
cat("done.\n")
# Show a preview of the knockoff stats
cat("Knockoff statistics:\n")
Stats

#######################################
## Save the knockoff test statistics ##
#######################################
out.file <- sprintf("%s/simulations_small/lasso/fold_%s/%s/stats_%s_%s.txt",
                    scratch, fold, resolution, experiment.name, experiment.numb)
Stats %>% write_delim(out.file, delim=" ")
cat(sprintf("Test statistics written to: %s\n", out.file))

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
Selections <- Stats %>% knockoff.filter(fdr=0.1, offset=1)

Selections <- Selections %>%
    select(Resolution, CHR, Group, W) %>% inner_join(Variants, by = c("CHR", "Group")) %>%
    group_by(Resolution, CHR, Group) %>%
    summarise(W=mean(W), BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n()/2,
              Causal=any(Causal)) %>%
    ungroup() %>%
    mutate(Method="Knockoffs")

# Give a preview of the selections
cat("Knockoff selections:\n")
Selections %>% print()

cat("Summary of knockoff selections:\n")
Selections %>% summarise(True=sum(Causal), False=sum(!Causal), FDP=False/(False+True)) %>% print()
