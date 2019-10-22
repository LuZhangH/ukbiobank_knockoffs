#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")
library(tidyverse)
library(kableExtra)

phenotype.name <- "platelet"

######################
## Load annotations ##
######################

scratch <- "/scratch/groups/candes"

######################
## Load discoveries ##
######################

# Load knockoffs discovereries
resolution.list <- c("Radj2", "Radj5", "Radj10", "Radj20", "Radj50", "Radj75", "Radj100")
phenotype.list <- c(phenotype.name)
Params <- expand.grid(Resolution=resolution.list, Phenotype=phenotype.list) %>% as_tibble()
Discoveries <- lapply(1:nrow(Params), function(idx) {
    resolution <- Params$Resolution[idx]
    phenotype <- Params$Phenotype[idx]
    knockoffs.file <- sprintf("%s/ukbiobank_tmp/discoveries/%s_knockoffs_%s.txt", scratch, phenotype, resolution)
    if(file.exists(knockoffs.file)) {
        Discoveries.knockoffs <- read_delim(knockoffs.file, delim=" ", col_types=cols()) %>%
            mutate(Phenotype=phenotype, Method="Knockoffs", Importance=W, Resolution=resolution) %>%
            dplyr::select(-c("W"))
    } else {
        Discoveries.knockoffs <- tibble()
    }
    return(Discoveries.knockoffs)
})
Discoveries <- do.call("rbind", Discoveries) %>% mutate(Resolution=as.character(Resolution))

Discoveries.hr <- Discoveries %>% filter(Resolution=="Radj100")

################################
## Cross reference with dbSNP ##
################################
library(AnnotationHub)
library(VariantAnnotation)

file.gz <- sprintf("%s/ukbiobank/annotations/00-All.vcf.gz", scratch)
tab <- TabixFile(file.gz)

## chr <- "12"
## pos <- 111885310
## snp <- "rs72650673"

n.read <- nrow(Discoveries.hr)
chr <- Discoveries.hr$CHR[1:n.read]
pos <- Discoveries.hr$BP.lead[1:n.read]
snp <- Discoveries.hr$SNP.lead[1:n.read]

df.original <- Discoveries.hr[1:n.read,] %>%
    transmute(CHR=as.integer(CHR), SNP=SNP.lead, RSPOS=as.integer(BP.lead))
df.original

rng <- GRanges(seqnames=chr, ranges=IRanges(start=pos, end=pos))
vcf_rng <- readVcf(tab, "hg19", param=rng)

# Look up annotations
snp.names <- rownames(info(vcf_rng))
df.annotations <- cbind(tibble(SNP=snp.names),info(vcf_rng)) %>% as_tibble()
df.annotations

# Cross-reference
df.annotations <- df.original %>% inner_join(df.annotations)
df.annotations

df.annotations %>% dplyr::select(-CHR,-SNP,-RSPOS,-GENEINFO) %>% summarise_all(sum)

# Genes
df.annotations %>% dplyr::select(CHR,SNP,GENEINFO)



# Look up annotations
snp.names <- rownames(info(vcf_rng))
df.annotations <- cbind(tibble(SNP=snp.names),info(vcf_rng)[,c(2,11,13,31,42)]) %>% as_tibble()
df.annotations <- df.original %>% inner_join(df.annotations)
df.annotations

df.annotations %>% dplyr::select(-CHR,-SNP,-RSPOS) %>% summarise_all(sum)

df.annotations %>% dplyr::select(-CHR,-RSPOS) %>%
    group_by(SNP) %>% mutate(referenced=any(PM,PMC,MUT,OM)) %>%
    ungroup() %>% summarise(referenced=sum(referenced))

head(rowRanges(vcf_rng), 3)
info(header(vcf_rng)) %>% as_tibble %>% print(n=100)

info(header(vcf_rng))$Description[c(11)]


# Gene and other info
info(vcf_rng)[,c(1,2,5,11,17)]

info(vcf_rng) %>% as_tibble()
