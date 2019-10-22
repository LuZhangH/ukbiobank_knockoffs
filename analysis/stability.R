#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

## Load packages
library(tidyverse)

## Other parameters
scratch <- "/scratch/PI/candes/ukbiobank_tmp"

## Default arguments
phenotype.list <- c("height","bmi","sbp","platelet","cvd","diabetes","hypothyroidism","respiratory","glaucoma")
resolution <- "Radj2"

Discoveries <- lapply(phenotype.list, function(phenotype) {
  ## Load original list of discoveries
  discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_%s.txt", scratch, phenotype, resolution)
  Discoveries.old <- read_delim(discoveries.file, delim=" ", col_types=cols()) %>% mutate(Seed=0)
  ## Load new list of discoveries
  seed <- 1
  discoveries.file <- sprintf("%s/discoveries/%s_knockoffs_%s_s%s.txt", scratch, phenotype, resolution, seed)
  Discoveries.new <- read_delim(discoveries.file, delim=" ", col_types=cols()) %>% mutate(Seed=1)
  ## Return results
  Discoveries.pheno <- rbind(Discoveries.old, Discoveries.new) %>% mutate(Phenotype=phenotype)
  return(Discoveries.pheno)
})
Discoveries <- do.call("rbind", Discoveries)

## Cross-reference
Discoveries.both <- inner_join(select(Discoveries.old,-Seed, -W), select(Discoveries.new,-Seed, -W))
## 
cat(sprintf("Replicated discoveries: %d out of %d (%.1f%%),\n",
            nrow(Discoveries.both), nrow(Discoveries.old), 100*nrow(Discoveries.both)/nrow(Discoveries.old)))
