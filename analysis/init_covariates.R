.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")
suppressMessages(library(tidyverse))
suppressMessages(library(superheat))
suppressMessages(library(bigsnpr))

scratch <- "/scratch/PI/candes/ukbiobank_tmp"

###############################
## Load the phenotype matrix ##
###############################

# Load header for phenotype table
pheno.header.file <- sprintf("%s/phenotypes/phenotypes_header.tab", scratch)
Phenotypes.header <- read_tsv(pheno.header.file, col_types=cols()) %>% mutate(Type=substring(Class,1,1))
col_types <- paste(Phenotypes.header$Type,collapse="")

# Load phenotype table
pheno.file <- sprintf("%s/phenotypes/phenotypes.tab", scratch)
Phenotypes <- read_tsv(pheno.file, col_types=col_types)

#################################################################
## Make list of response variables for analysis and covariates ##
#################################################################

phenotypes.exclude <- c("FID", "IID", "sex", "age", "age.sq", "PC.1", "PC.2", "PC.3", "PC.4", "PC.5")

# List of default covariates
covariates.default <- c("sex", "age", "age.sq")

# List of custom covariates
covariates.custom <- tribble(
  ~Name, ~Covariates,
  "height",   c("sex", "age"),
)

# Make list of response variables and corresponding types
Analysis <- Phenotypes.header %>% select(-Type) %>%
    filter(! Name %in% phenotypes.exclude)
# Assign list of covariates
Analysis$Covariates <- list(covariates.default)
# Modify custom covariates
for(i in 1:nrow(covariates.custom)) {
    name <- covariates.custom$Name[i]
    Analysis$Covariates[which(Analysis$Name==name)] <- covariates.custom$Covariates[i]
}

cat("List of response variables for analysis:\n")
Analysis %>% rowwise() %>% mutate(Covariates=paste( unlist(Covariates), collapse=',')) %>% print()

analysis.file <- sprintf("%s/phenotypes/analysis.tab", scratch)
Analysis %>%
    rowwise() %>% mutate(Covariates=paste( unlist(Covariates), collapse=',')) %>% 
    write_tsv(analysis.file)
cat(sprintf("Writing list phenotypes and covariates for analysis in:\n %s\n", analysis.file))
    
# Write lists on files for BOLT-LMM
covariates.dir <- sprintf("%s/analysis/bolt", scratch)

for(i in 1:nrow(Analysis)) {
    pheno.name <- Analysis$Name[i]
    covariates.file <- sprintf("%s/%s_covariates.txt", covariates.dir, pheno.name)
    
    cat(sprintf("Writing list of BOLT covariates for phenotype %s in:\n %s\n", pheno.name, covariates.file))
    sink(covariates.file)
    
    for(covariate.name in Analysis$Covariates[[i]]) {
        pheno.row <- which(Phenotypes.header$Name==covariate.name)
        covariate.type <- Phenotypes.header$Type[pheno.row]
        if(covariate.type %in% c("i","n")) {
            cat(sprintf("--qCovarCol=%s \n", covariate.name, covariate.type))
        } else {
            cat(sprintf("--covarCol=%s \n", covariate.name, covariate.type))
        }
    }
    sink()
}
