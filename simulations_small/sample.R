#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

library(tidyverse)
library(snpStats)
library(bigsnpr)

scratch <- "/scratch/groups/candes/ukbiobank_tmp"

set.seed(2019)

#########################
## List of experiments ##
#########################
n.folds <- 1

# Define table of possible architectures
#name.list <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "sevenps")
#n.loci.list <- c(rep(100,3), rep(500,3), rep(1000,4), 500)
#s.loci.list <- c(rep(c(1,2,5),3), 10, 5)

name.list <- c("sevenps")
n.loci.list <- c(500)
s.loci.list <- c(1)

Architecture <- tibble(Name=name.list, Locus.number=n.loci.list, Locus.size=s.loci.list) %>%
    mutate(Causal=Locus.number*Locus.size) %>%
    select(Name, Causal, everything()) %>%
    arrange(Locus.number, Locus.size)

# Save architecture table
table.file <- sprintf("%s/simulations_small/phenotypes/architecture.txt", scratch)
Architecture %>% write_delim(table.file, delim=" ")

###########################################
## Load list of variants and individuals ##
###########################################

# Load fam file
fam.file <- sprintf("%s/knockoffs/Radj10_K50/ukb_gen_chr%d.fam", scratch, 22)
Subjects <- read_delim(fam.file, delim=" ", col_types=cols(),
                       col_names = c("FID", "IID", "X1", "X2", "Sex", "X3"))

# Load list of original variants
Variants <- tibble()
for(chr in seq(1,22)) {
    cat(sprintf("Loading list of variants on chromosome %d... ", chr))
    key.file <- sprintf("%s/knockoffs/Radj10_K50/ukb_gen_chr%d.key", scratch, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% filter(Knockoff==FALSE) %>%
        select(Chr, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants <- rbind(Variants, Variants.chr)
    cat("done.\n")
}

################################
## Define signal architecture ##
################################

choose_causal <- function(n.loci, n.signals) {

    # Check that the signals can be split exactly into the desired number of loci
    locus.size <- n.signals/n.loci
    if(locus.size!=round(locus.size)) {
        cat(sprintf("Error: %d is not divisible by %d.\n", n.signals, n.loci))
        return(NULL)
    }

    # Allocate causal SNPs to chromosomes
    Variants.num <- Variants %>% group_by(CHR) %>% summarise(N=n())
    Variants.num <- Variants.num %>% mutate(Share = N/sum(Variants.num$N), Loci = round(Share*n.loci))
    # Add signals if rounding errors occurred
    if(sum(Variants.num$Loci) != n.loci) {
        Variants.num$Loci[1] = Variants.num$Loci[1]+(n.loci-sum(Variants.num$Loci))
    }
    n.loci.chr <- Variants.num$Loci
    cat(sprintf("Divided %d causal loci into %d chromosomes.\n", sum(n.loci.chr),nrow(Variants.num)))

    # Distance between causal SNPs (in kB)
    causal.dist <- round(sort(runif(locus.size, 0, 100)),2)

    # Pick causal SNPs
    Variants.causal <- tibble()
    for(chr in seq(1,22)) {
        cat(sprintf("Assembling %d clumps of size %d for chromosome %d...\n", n.loci.chr[chr], locus.size, chr))
        Variants.chr <- Variants %>% filter(CHR==chr)
        causal.idx <- round(seq(100,round(0.95*nrow(Variants.chr)), length.out=n.loci.chr[chr]))
        bp.first <- Variants.chr$BP[causal.idx]
        bp.next <- as.numeric(sapply(causal.dist[-1], function(dkb) {bp.first + 1000 * dkb}))
        bp.next <- matrix(bp.next, ncol=length(causal.dist)-1)
        locus.idx <- seq(length(causal.idx))
        for(j in 1:length(bp.first)) {
            if(length(causal.dist)>1) {
                for(k in 1:ncol(bp.next)) {
                    candidates <- order(abs(Variants.chr$BP-bp.next[j,k]))
                    candidates <- candidates[which(candidates>causal.idx[j])]
                    candidates <- setdiff(candidates, causal.idx)
                    causal.idx <- c(causal.idx, candidates[1])
                    locus.idx <- c(locus.idx, j)
                }
            }
        }
        Variants.causal.chr <- Variants.chr[causal.idx,]
        Variants.causal.chr$Locus <- locus.idx
        Variants.causal.chr <- Variants.causal.chr %>% arrange(BP) %>% distinct()
        Variants.causal <- rbind(Variants.causal, Variants.causal.chr)
    }

    # Assign effect sizes
    Variants.causal <- Variants.causal %>%
        mutate(Locus=rep(seq(1,n.loci),each=locus.size)) %>%
        mutate(Sign=rep(2*rbinom(nrow(Variants.causal)/locus.size,1,0.5)-1,each=locus.size)) %>%
        mutate(Scale=runif(nrow(Variants.causal),0.1,1.9))

    return(Variants.causal)
}

Variants.causal.all <- lapply(1:nrow(Architecture), function(exp.num) {
    # Choose causal variants for this experiment
    exp.name <- Architecture$Name[exp.num]
    n.loci <- Architecture$Locus.number[exp.num]
    n.signals <- Architecture$Causal[exp.num]
    Variants.causal <- choose_causal(n.loci, n.signals) %>% mutate(Experiment=exp.name)
    # Save list of causal variants
    causal.file <- sprintf("%s/simulations_small/phenotypes/%s_causal.txt", scratch, exp.name)
    write_delim(Variants.causal, causal.file, delim=" ")
    cat(sprintf("Written list of %d causal variants for experiment %s on:\n %s\n",
                nrow(Variants.causal), exp.name, causal.file))
    return(Variants.causal)
})
Variants.causal.all <- do.call("rbind", Variants.causal.all)

########################
## Load the genotypes ##
########################

rds.file <- sprintf("%s/augmented_data_big/ukb_gen_Radj10.rds", scratch)
if(file.exists(rds.file)){
    cat(sprintf("Found FBM in %s.\n", rds.file))
} else {
    cat(sprintf("Could not find FBM in %s.\n", rds.file))
    quit()
}
# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
obj.bigSNP <- snp_attach(rds.file)
cat("done.\n")
# Compute scaling factor for the genotypes
scale.file <- sprintf("%s/augmented_data_big/ukb_gen_Radj10_scale.txt", scratch)
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

##########################
## Load the confounders ##
##########################

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
# Assign relative effects sizes to each principal component
beta.Z <- 1/seq(ncol(Z)) * rep(c(+1,-1), length.out=ncol(Z))

##################
## Sample Y | X ##
##################

sample_phenotypes_real <- function(obj.bigSNP, exp.name, Variants.causal.all,
                                   h2.values, h2.values.confounding=NULL,
                                   Architecture) {
    Variants.causal <- Variants.causal.all %>% filter(Experiment==exp.name)
    
    # Get aliases for useful slots
    G <- obj.bigSNP$genotypes
    CHR <- obj.bigSNP$map$chromosome
    POS <- obj.bigSNP$map$physical.pos
    SNP <- obj.bigSNP$map$marker.ID

    # Extract causal columns and standardize them
    col.idx <- which(SNP %in% Variants.causal$SNP)
    cat(sprintf("Scaling and centering genotypes for %d causal variants... ", length(col.idx)))
    X <- G[,col.idx]
    X <- sweep(X, 2, colMeans(X), FUN = '-')
    X <- sweep(X, 2, as.numeric(scaling.factors[col.idx]), FUN = '/')
    cat("done.\n")
    
    # Compute phenotype mean vector
    cat(sprintf("Computing phenotype means... "))
    beta <- Variants.causal$Sign * Variants.causal$Scale / sqrt(nrow(X))
    if(is.null(h2.values.confounding)) {
        # No confounding 
        y.mean <- X %*% beta

        # Estimate heritability and choose interesting signal amplitudes
        cat(sprintf("Computing signal amplitudes... "))
        y.mean.var <- var(y.mean)
        amplitude.curve <- seq(0.1,100, length.out=10000)
        h2.curve <- sapply(amplitude.curve, function(amplitude) {
            h2 <- amplitude^2*y.mean.var / (amplitude^2*y.mean.var + 1)
            return(h2)
        })
        amplitude.values <- sapply(h2.values, function(h2) {
            amplitude <- amplitude.curve[which.min(abs(h2.curve-h2))]
            return(round(amplitude,3))
        })
        cat("done.\n")
        
        # Sample Y | X
        Phenotypes <- obj.bigSNP$fam %>% as_tibble() %>% transmute(FID=family.ID, IID=sample.ID)
        for(a.idx in seq(length(amplitude.values))) {
            for(fold in 1:n.folds) {
                amplitude <- amplitude.values[a.idx]
                cat(sprintf("Sampling Y|X for amplitude %.3f (fold %d).\n", amplitude, fold))
                y <- round(amplitude * y.mean + rnorm(nrow(X)),4)
                phenotype.colnames <- colnames(Phenotypes)
                y.name <- sprintf("Y_a%s_%s", a.idx, fold)
                phenotype.colnames <- c(phenotype.colnames, y.name)
                Phenotypes <- Phenotypes %>% mutate(Y = y)
                colnames(Phenotypes) <- phenotype.colnames
            }
        }
    } else {         
        # Include confounding 
        h2.values.total <- h2.values.confounding + h2.values[1]
        
        # Choose amplitude for sparse signal
        y.mean.sparse <- X %*% beta
        cat(sprintf("Computing signal amplitudes for sparse signals... "))
        y.mean.sparse.var <- var(y.mean.sparse)
        amplitude.curve <- seq(0.1,100, length.out=10000)
        h2.curve <- sapply(amplitude.curve, function(amplitude) {
            h2 <- amplitude^2*y.mean.sparse.var / (amplitude^2*y.mean.sparse.var + 1)
            return(h2)
        })
        amplitude.values <- sapply(h2.values, function(h2) {
            amplitude <- amplitude.curve[which.min(abs(h2.curve-h2))]
            return(round(amplitude,4))
        })
        cat("done.\n")

        # Choose amplitude for confounders
        cat(sprintf("Computing signal amplitudes for confounders... "))
        amplitude.curve <- seq(0,10, length.out=10000)
        var.xb <- var(X %*% beta) %>% as.numeric()
        var.zb <- var(Z %*% beta.Z) %>% as.numeric()
        cov.xz <- cov(X %*% beta, Z %*% beta.Z) %>% as.numeric()
        h2.curve <- sapply(amplitude.curve, function(amplitude) {
            y.mean.var <- amplitude.values[1]^2*var.xb + amplitude^2*var.zb
            y.mean.var <- y.mean.var +  2*amplitude.values[1]*amplitude * var.xb*var.zb
            h2 <- y.mean.var / (y.mean.var + 1)
            return(h2)
        })
        amplitude.values.Z <- sapply(h2.values.total, function(h2) {
            amplitude <- amplitude.curve[which.min(abs(h2.curve-h2))]
            return(round(amplitude,4))
        })
        cat("done.\n")
        
        # Sample Y | X
        Phenotypes <- obj.bigSNP$fam %>% as_tibble() %>% transmute(FID=family.ID, IID=sample.ID)
        for(a.idx in seq(length(amplitude.values.Z))) {
            y.mean <- amplitude.values[1] * X %*% beta + amplitude.values.Z[a.idx] * Z %*% beta.Z
            for(fold in 1:n.folds) {
                amplitude.Z <- amplitude.values.Z[a.idx]
                cat(sprintf("Sampling Y|X for amplitude.Z %.3f (fold %d).\n", amplitude.Z, fold))
                y <- round(y.mean + rnorm(nrow(X)),4)
                phenotype.colnames <- colnames(Phenotypes)
                y.name <- sprintf("Y_a%s_%s", a.idx, fold)
                phenotype.colnames <- c(phenotype.colnames, y.name)
                Phenotypes <- Phenotypes %>% mutate(Y = y)
                colnames(Phenotypes) <- phenotype.colnames
            }
        }
        # Store the H2 and amplitude values for Z
        h2.values <- h2.values.confounding
        amplitude.values <- amplitude.values.Z
    }
    
    cat("done.\n")

    # Make sure that the rows of the fam file match the rows of the phenotypes
    Phenotypes <- Phenotypes %>% right_join(select(Subjects,FID,IID), by = c("FID", "IID"))

    # Save phenotypes
    pheno.file <- sprintf("%s/simulations_small/phenotypes/%s_phenotypes.tab", scratch, exp.name)
    write_delim(Phenotypes, pheno.file, delim="\t")
    cat(sprintf("Written phenotypes on:\n %s\n", pheno.file))

    # Update architecture
    Architecture.new <- tibble(Name=exp.name, H2=round(h2.values,3),
                               Amplitude=seq(length(amplitude.values)), Amplitude.num=amplitude.values)        
    return(Architecture.new)
}

liability.model <- function(X,beta,case.fraction) {
    L <- X %*% beta + rnorm(nrow(X))    
    L.threshold <- quantile(L, 1-case.fraction/(1+case.fraction))
    y <- 1+as.numeric(L>L.threshold)
    return(y)
}

sample_phenotypes_binary <- function(obj.bigSNP, exp.name, Variants.causal.all,
                                     h2.values, h2.values.confounding=NULL, Architecture) {
    Variants.causal <- Variants.causal.all %>% filter(Experiment==exp.name)
    
    # Get aliases for useful slots
    G <- obj.bigSNP$genotypes
    CHR <- obj.bigSNP$map$chromosome
    POS <- obj.bigSNP$map$physical.pos
    SNP <- obj.bigSNP$map$marker.ID

    # Extract causal columns and standardize them
    col.idx <- which(SNP %in% Variants.causal$SNP)
    cat(sprintf("Scaling and centering genotypes for %d causal variants... ", length(col.idx)))
    X <- G[,col.idx]
    X <- sweep(X, 2, colMeans(X), FUN = '-')
    X <- sweep(X, 2, as.numeric(scaling.factors[col.idx]), FUN = '/')
    cat("done.\n")
    
    # Compute phenotype mean vector
    cat(sprintf("Computing phenotype means... "))
    case.fraction <- 0.0015
    beta <- Variants.causal$Sign*Variants.causal$Scale / sqrt(nrow(X))
    
    # Estimate heritability and choose interesting signal amplitudes
    cat(sprintf("Computing signal amplitudes... "))
    L.mean <- X %*% beta
    L.mean.var <- var(L.mean)
    amplitude.curve <- seq(0.1,100, length.out=10000)
    h2.curve <- sapply(amplitude.curve, function(amplitude) {
        h2 <- amplitude^2*L.mean.var / (amplitude^2*L.mean.var + 1)
        return(h2)
    })
    amplitude.values <- sapply(h2.values, function(h2) {
        amplitude <- amplitude.curve[which.min(abs(h2.curve-h2))]
        return(round(amplitude,3))
    })
    cat("done.\n")

    # Sample Y | X
    Phenotypes <- obj.bigSNP$fam %>% as_tibble() %>% transmute(FID=family.ID, IID=sample.ID)
    for(a.idx in seq(length(amplitude.values))) {
        for(fold in 1:n.folds) {
            amplitude <- amplitude.values[a.idx]
            cat(sprintf("Sampling Y|X for amplitude %.3f (fold %d).\n", amplitude, fold))
            y <- liability.model(X, amplitude*beta, case.fraction)
            phenotype.colnames <- colnames(Phenotypes)
            y.name <- sprintf("Y_a%s_%s", a.idx, fold)
            phenotype.colnames <- c(phenotype.colnames, y.name)
            Phenotypes <- Phenotypes %>% mutate(Y = y)
            colnames(Phenotypes) <- phenotype.colnames
        }
    }

    # Make sure that the rows of the fam file match the rows of the phenotypes
    Phenotypes <- Phenotypes %>% right_join(select(Subjects,FID,IID), by = c("FID", "IID"))

    # Save phenotypes
    pheno.file <- sprintf("%s/simulations_small/phenotypes/%s_phenotypes.tab", scratch, exp.name)
    write_delim(Phenotypes, pheno.file, delim="\t")
    cat(sprintf("Written phenotypes on:\n %s\n", pheno.file))

    # Update architecture
    Architecture.new <- tibble(Name=exp.name, H2=round(h2.values,3),
                               Amplitude=seq(length(amplitude.values)), Amplitude.num=amplitude.values)        
    return(Architecture.new)
}

Architecture.full <- tibble()
for(exp.name in name.list) {

    if(exp.name == "sevenps") {
        h2.values <- c(0.35)
        h2.values.confounding <- c(0, 0.01, 0.05, 0.075, round(seq(0.1, 0.35, length.out=6),3))
    } else {
        h2.values <- c(0.01, 0.05, 0.075, seq(0.1, 0.7, by=0.1))
        h2.values.confounding <- NULL
    }

    if(exp.name=="zero") {
        Architecture.new <- sample_phenotypes_binary(obj.bigSNP, exp.name, Variants.causal.all,
                                                     h2.values,
                                                     h2.values.confounding=h2.values.confounding,
                                                     Architecture)
    } else {
        Architecture.new <- sample_phenotypes_real(obj.bigSNP, exp.name, Variants.causal.all,
                                                   h2.values,
                                                   h2.values.confounding=h2.values.confounding,
                                                   Architecture)
    }

    if(nrow(Architecture.full)>0) {
        Architecture.new <- inner_join(Architecture, Architecture.new, by = "Name")
        Architecture.full <- rbind(Architecture.full, Architecture.new) %>% as_tibble()
    } else {
        Architecture.full <- inner_join(Architecture, Architecture.new, by = "Name")
    }
    Architecture.full %>%
        mutate(Name=factor(Name, levels=name.list, labels=name.list)) %>%
        arrange(Name, Locus.number, Locus.size, H2) %>% select(Name, Amplitude, everything()) %>%
        write_delim(table.file, delim=" ")
    cat(sprintf("Updated architecture file:\n %s\n",table.file))
}
