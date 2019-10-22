#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Input arguments
chr <- 22
seed <- 123
mask <- "01"
    
# Load libraries
suppressMessages(library(snpStats))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(SNPknock))
source("../../utils/util.R")
options(scipen = 999)

# Location of data
tmp.dir <- "/scratch/PI/candes/ukbiobank_tmp/fastphase"

load.results <- function(n, K.list){
    fp.dir <- paste(tmp.dir, "/CV_", n, sep="")
    # Load haplotypes from INP file
    basename <- paste(fp.dir, "/input/ukb_hap_chr", chr, "_test", sep="")
    Haplotypes <- read.inp(basename, progress=TRUE) %>% as_tibble()
    #H <- as.matrix(dplyr::select(Haplotypes, -Subject))
    H <- as.matrix(Haplotypes)

    # Load masked haplotypes from INP file
    basename <- paste(fp.dir, "/input/mask",mask,"/ukb_hap_chr", chr, "_test", sep="")
    Haplotypes.masked <- read.inp(basename, progress=TRUE)
    H.masked <- as.matrix(Haplotypes.masked)

    # Load imputed haplotypes from INP file
    Results <- sapply(K.list, function(K) {
        # Load haplotypes imputed with this value of K
        basename <- paste(fp.dir, "/K",K,"/mask",mask,"/ukb_hap_chr", chr, "_imputed", sep="")
        if(!file.exists(paste(basename,".inp",sep=""))) {
            return(c(n,K,NA))
        }
        Haplotypes.imputed <- read.inp(basename, progress=TRUE)
        H.imputed <- as.matrix(Haplotypes.imputed)
        # Compute imputation error
        error <- mean(H.imputed[is.na(H.masked)]!=H[is.na(H.masked)])
        return(c(n,K,error))
    })
    Results <- as.tibble(t(Results))
    colnames(Results) <- c("n", "K", "error")

    return(Results)
}

# Load results from simulations
K.list <- c(1,2,5,10,20,30,50,75,100,150,200)
Results <- load.results(100, K.list)
for( n in c(300,1000,3000,10000,30000,100000)) {
    cat(sprintf("Loading results with n = %d\n", n))
    Results <- rbind(Results,load.results(n, K.list))
}

# Load results from full data
for(K in c(1,2,5,10,15,20,30,50,75,100)) {
    filename <- paste(tmp.dir, "/test_imputed/phased_K", K, "/ukb_hap_chr", chr, "_error.txt", sep="")
    res <- suppressMessages(read_delim(filename, delim="\t"))
    res <- res %>% dplyr::select(c("n", "K", "error")) %>% mutate(n=as.integer(349119))
    Results <- rbind(Results,res)
}

# Save results
out.dir <- "/scratch/PI/candes/ukbiobank_tmp/reports"
out.file <- sprintf("%s/goodnessofit.txt", out.dir)
Results %>% write_delim(out.file, delim=" ")
cat(sprintf("Results saved in:\n %s\n", out.file))

# Plot results
Results %>% mutate(n=as.ordered(n)) %>% ggplot(aes(x=K, y=error, color=n)) +
    geom_point() + geom_line() + scale_colour_hue(h = c(100, 350)) + theme_bw() + ylim(0,0.12)
