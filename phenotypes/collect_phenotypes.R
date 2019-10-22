library(tidyverse)

# Data storage
pheno.dir <- "/scratch/PI/candes/ukbiobank/phenotypes"
pheno.raw <- sprintf("%s/ukb25261.csv", pheno.dir)
scratch <- "/scratch/PI/candes/ukbiobank_tmp"

#############################
## Specify data to extract ##
#############################
# Parameters
n.PCs <- 5

# List of quantitative measures to extract
Measures <- tribble(
  ~Name, ~Type, ~Code,
  "sex",           "categorical",  "31-0.0",
  "age",           "continuous",   "21022-0.0",
  "height",        "continuous",   "50-0.0",
  "bmi",           "continuous",   "21001-0.0",
  "waist",         "continuous",   "48-0.0",
  "waist",         "continuous",   "49-0.0",
  "dbp",           "continuous",   "4079-0.0",
  "dbp",           "continuous",   "4079-0.1",
  "sbp",           "continuous",   "4080-0.0",
  "sbp",           "continuous",   "4080-0.1",
  "platelet",      "continuous",   "30080-0.0",
  "redcell",       "continuous",   "30010-0.0",
  "whitecell",     "continuous",   "30000-0.0",
  "eosinophil",    "continuous",   "30150-0.0",
  "redwidth",      "continuous",   "30070-0.0",
  "plateletwidth", "continuous",   "30110-0.0",
)

# List of diseases to extract
Diseases <- tribble(
  ~Name, ~Code, ~Prefix,
    "cvd", c(1065,1066,1067,1068,1081,1082,1083,1425,1473,1493), "20002-0.",
    "diabetes", 1220, "20002-0.",
    "hypothyroidism", 1226, "20002-0.",
    "respiratory", c(1111,1112,1113,1114,1115,1117,1413,1414,1415,1594), "20002-0.",
)
Diseases <- Diseases %>% mutate(Type="categorical")
                             
# List of diagnoses to extract
Diagnoses <- tribble(
  ~Name, ~Code, ~Prefix,
    "glaucoma", c("H400","H401","H402","H403","H404","H405","H406","H408","H409",
                  "H42", "H420", "H428"), "41202-0.",
)
Diagnoses <- Diagnoses %>% mutate(Type="categorical")

# Columns containing principal components
PC.cols <- paste("22009-0.", seq(1,n.PCs), sep="")
# Columns containing various measurements
mea.cols <- Measures$Code
# Columns containing self-reported disease status
dis.cols <- paste("20002-0.", seq(0,32), sep="")
# Columns containing diagnoses (not used)
dia.cols <- paste("41202-0.", seq(0,379), sep="")

# Form complete list of columns
cols.extract <- unique(c("eid", PC.cols, mea.cols, dis.cols, dia.cols))
cat(sprintf("The following columns will be extracted from the original phenotype file:\n"))
print(cols.extract)

##########################################################
## Load raw data and initialize refined phenotype table ##
##########################################################

# Read complete file (this will take a couple of minutes)
Data.raw <- read_csv(pheno.raw, guess_max=500000)

# Check whether all requested phenotypes are available from the original file
not.found <- which(! cols.extract %in% colnames(Data.raw))
if(length(not.found) > 0) {
    cat(sprintf("%d out of %d columns were not found in the original phenotype file:\n",
                length(not.found), length(cols.extract)))
    print(cols.extract[not.found])
} else {
    cat(sprintf("All %d columns were found in the original phenotype file.\n", length(cols.extract)))
}

# Extract requested columns
Data <- Data.raw %>% select(cols.extract)

# Initialize phenotype table
Phenotypes <- select(Data, "eid")

###############################
## Extract physical measures ##
###############################

# Extract columns corresponding to physical measures
Data.measures <- Data %>% select("eid", Measures$Code)
colnames(Data.measures) <- c("eid", Measures$Name)

# Combine columns with the same name
Data.measures <- sapply(colnames(Data.measures), function(colname) {
    col.idx <- which(colnames(Data.measures)==colname)
    col.mean <- rowMeans(Data.measures[,col.idx], na.rm=T)
    col.mean[is.nan(col.mean)] <- NA
    return(col.mean)
}) %>% as_tibble() %>% select(colnames(Data.measures))

# Sex encoding (0 = missing, 1 = male, 2 = female)
Data.measures <- Data.measures %>% mutate(sex=sex+1)
sex.is.male <- which(Data.measures$sex==1)
sex.is.female <- which(Data.measures$sex==0)
sex.is.missing <- which(is.na(Data.measures$sex))
Data.measures$sex[sex.is.male] <- 1
Data.measures$sex[sex.is.female] <- 2
Data.measures$sex[sex.is.missing] <- 0

# Insert these columns into the new phenotype table
Phenotypes <- left_join(Phenotypes, Data.measures, by=c("eid"))

######################
## Extract diseases ##
######################

# Extract binary self-reported diseases
for(i in 1:nrow(Diseases)) {
    disease.name <- Diseases$Name[[i]]
    disease.code <- Diseases$Code[[i]]
    disease.prefix <- Diseases$Prefix[[i]]
    cat(sprintf("Collecting phenotype '%s' (%s) from fields '%s' ...\n",
                disease.name, disease.code, disease.prefix))
    matrix.disease <- Data %>% select(starts_with(disease.prefix)) %>%
        map( function(x) x %in% disease.code ) %>% as_tibble() %>% as.matrix
    Data.disease <- tibble(eid=Data$eid, disease=1+as.logical(rowSums(matrix.disease)))
    Phenotypes[disease.name] <- Data.disease$disease
    disease.table <- table(Data.disease$disease)
    cat(sprintf("Found %d cases and %d controls.\n", disease.table[2], disease.table[1]))
}

#######################
## Extract diagnoses ##
#######################

# Extract binary diagnosed diseases
for(i in 1:nrow(Diagnoses)) {
    diagnosis.name <- Diagnoses$Name[[i]]
    diagnosis.code <- Diagnoses$Code[[i]]
    diagnosis.prefix <- Diagnoses$Prefix[[i]]
    cat(sprintf("Collecting phenotype '%s' (%s) from fields '%s' ...\n",
                diagnosis.name, diagnosis.code, diagnosis.prefix))
    matrix.diagnosis <- Data %>% select(starts_with(diagnosis.prefix)) %>%
        map( function(x) x %in% diagnosis.code ) %>% as_tibble() %>% as.matrix
    Data.diagnosis <- tibble(eid=Data$eid, diagnosis=1+as.logical(rowSums(matrix.diagnosis)))
    Phenotypes[diagnosis.name] <- Data.diagnosis$diagnosis
}

##################################
## Extract principal components ##
##################################

# Extract PCs
Data.PC <- Data %>% select(eid, PC.cols)
names(Data.PC) <- gsub(x = names(Data.PC), pattern = "22009-0.", replacement = "PC.")

# Insert these columns into the new phenotype table
Phenotypes <- left_join(Phenotypes, Data.PC, by=c("eid"))

##############################
## Summarise and save table ##
##############################

# Rename sample id column
Phenotypes <- Phenotypes %>% mutate(FID=as.integer(eid), IID=as.integer(eid)) %>%
    select(FID, IID, everything(), -eid)

# Include transformed variables
Phenotypes <- Phenotypes %>% mutate(age.sq = age^2)

# Convert column types to the correct format
for(i in 1:nrow(Measures)) {
    col.name <- Measures$Name[i]
    col.type <- Measures$Type[i]
    if(col.type=="continuous") {
        Phenotypes[[col.name]] <- as.numeric(Phenotypes[[col.name]])
    } else {
        Phenotypes[[col.name]] <- as.factor(Phenotypes[[col.name]])
        cat(sprintf("Converted %s to factor type\n", col.name))
    }
}
for(i in 1:nrow(Diseases)) {
    col.name <- Diseases$Name[i]
    Phenotypes[[col.name]] <- as.factor(Phenotypes[[col.name]])
    cat(sprintf("Converted %s to factor type\n", col.name))
}
for(i in 1:nrow(Diagnoses)) {
    col.name <- Diagnoses$Name[i]
    Phenotypes[[col.name]] <- as.factor(Phenotypes[[col.name]])
    cat(sprintf("Converted %s to factor type\n", col.name))
}

# Count cases
cat(sprintf("Number of disease cases:\n"))
Phenotypes.binary <- select(Phenotypes, Diseases$Name)
Phenotypes.binary %>% summarise_all(function(x) {mean(x==2)})

# Save phenotype header
out.header.file <- sprintf("%s/phenotypes/phenotypes_header.tab", scratch)
Phenotypes.header <- tibble(Name=colnames(Phenotypes), Class=sapply(Phenotypes, class))
write_delim(Phenotypes.header, out.header.file, delim="\t")
cat(sprintf("Saved header of phenotype matrix in:\n %s\n", out.header.file))

# Save phenotypes
out.file <- sprintf("%s/phenotypes/phenotypes.tab", scratch)
write_delim(Phenotypes, out.file, delim="\t")
cat(sprintf("Saved phenotype matrix in:\n %s\n", out.file))

#############################################################
## Cross-reference with list of individuals that passed QC ##
#############################################################

# Load list of subjects that passed QC
qc.file <- sprintf("%s/QC_output/individuals_QC.txt", scratch)
Subjects <- read_tsv(qc.file, col_types=cols(), col_names = c("FID", "IID"))

Phenotypes.qc <- Phenotypes %>% inner_join(Subjects, by = c("FID", "IID"))
out.file <- sprintf("%s/phenotypes/phenotypes_qc.tab", scratch)
write_delim(Phenotypes.qc, out.file, delim="\t")
cat(sprintf("Saved phenotype matrix for individuals that passed QC in:\n %s\n", out.file))
