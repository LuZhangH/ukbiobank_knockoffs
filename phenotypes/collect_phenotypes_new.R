library(tidyverse)

# Data storage
pheno.dir <- "/scratch/PI/candes/ukbiobank/phenotypes"
pheno.raw <- sprintf("%s/ukb26877.csv", pheno.dir)
scratch <- "/scratch/PI/candes/ukbiobank_tmp"

###################
## Load raw data ##
###################

# Read complete file (this will take a few minutes)
Data.raw <- read_csv(pheno.raw, col_types=cols(.default=col_character()), progress=TRUE)

#############################
## Specify data to extract ##
#############################

# List of quantitative measures to extract
Measures <- tribble(
  ~Name, ~Type, ~Code,
  "glaucoma",      "continuous",  c("4689-0.0","4689-1.0","4689-2.0")
)

# Columns containing various measurements
mea.codes <- unlist(Measures$Code)
mea.names <- Measures$Name

# Form complete list of columns
cols.extract <- unique(c("eid", mea.codes))
cat(sprintf("The following columns will be extracted from the original phenotype file:\n"))
print(cols.extract)

#######################
## Process raw data  ##
#######################

Data <- Data.raw %>% select(eid, mea.codes) %>% type_convert()

# Process glaucoma
Data <- Data %>%
    mutate(`4689-0.0`=replace(`4689-0.0`, `4689-0.0`<0, NA)) %>%
    mutate(`4689-1.0`=replace(`4689-1.0`, `4689-1.0`<0, NA)) %>%
    mutate(`4689-2.0`=replace(`4689-2.0`, `4689-2.0`<0, NA)) %>%
    mutate(glaucoma=pmin(`4689-0.0`,`4689-1.0`,`4689-2.0`,na.rm=T)) %>%
    mutate(glaucoma=glaucoma > 30) %>% 
    mutate(glaucoma=replace(glaucoma, glaucoma==TRUE, 2)) %>%
    mutate(glaucoma=replace(glaucoma, glaucoma==FALSE, 1)) %>%
    mutate(glaucoma=replace(glaucoma, is.na(glaucoma), 1)) %>%
    mutate(glaucoma=as.integer(glaucoma)) %>%
    select(eid, glaucoma)

##############################
## Summarise and save table ##
##############################

# Rename sample id column
Phenotypes <- Data %>% mutate(FID=as.integer(eid), IID=as.integer(eid)) %>%
    select(FID, IID, everything(), -eid)

# Save phenotypes
out.file <- sprintf("%s/phenotypes/phenotypes_new.tab", scratch)
write_delim(Phenotypes, out.file, delim="\t")
cat(sprintf("Saved phenotype matrix in:\n %s\n", out.file))
