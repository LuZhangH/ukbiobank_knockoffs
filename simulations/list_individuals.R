#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

library(tidyverse)

scratch <- "/scratch/groups/candes/ukbiobank_tmp"

# Load list of all subjects
fam.file <- "/scratch/groups/candes/ukbiobank_tmp/simulations/subjects/ukb_gen.fam"
col_names <- c("FID", "IID", "X1", "X2", "X3", "X4")
FAM <- read_delim(fam.file, delim=" ", col_names=col_names)

# Load list of subjects that passed QC
subjects.file <- "/scratch/groups/candes/ukbiobank_tmp/augmented_data_big/ukb_gen_Radj100.fam"
Subjects <- read_delim(subjects.file, delim=" ", col_names=col_names)

# Make list of removed individuals
Remove <- anti_join(FAM, Subjects, by = c("FID", "IID", "X1", "X2", "X3", "X4"))
remove.file <- "/scratch/groups/candes/ukbiobank_tmp/simulations/subjects/ukb_gen_remove.txt"
Remove %>% write_delim(remove.file, delim=" ")
