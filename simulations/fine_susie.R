#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")
library(devtools)
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))
suppressMessages(library(bigsnpr))
source("utils/fine_mapping.R")
source("../utils/utils_clumping.R")

# Default arguments
experiment.name <- "six"
amplitude <- 10
fold <- 1
clumping <- "0.00000005"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
experiment.name <- as.character(args[1])
amplitude <- as.character(args[2])
fold <- as.character(args[3])
clumping <- as.character(args[4])

# Other parameters
scratch <- "/scratch/PI/candes/ukbiobank_tmp"
chr.list <- seq(1,22)
resolution <- "Radj100"
resolution.list <- c(resolution)
pheno.name <- sprintf("Y_a%s_%s", amplitude, fold)

###########################
## Load list of variants ##
###########################

# Load list of variants
Variants <- lapply(chr.list, function(chr) {
    key.file <- sprintf("%s/knockoffs/%s_K50/ukb_gen_chr%d.key", scratch, resolution, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants.chr$Resolution <- resolution
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants)

# Load list of causal variants
causal.file <- sprintf("%s/simulations/phenotypes/%s_causal.txt", scratch, experiment.name)
Causal <- read_delim(causal.file, delim=" ", col_types=cols()) %>% select(CHR, BP, Locus, Sign, Scale)
Variants <- Variants %>%
    left_join(mutate(Causal,Causal=TRUE),by = c("CHR", "BP")) %>%
    mutate(Causal=replace_na(Causal,FALSE)) %>%
    mutate(Causal=ifelse(Knockoff, FALSE, Causal))

# Extract list of original variants
Variants.original <- Variants %>%
    filter(Resolution==resolution.list[1], Knockoff==FALSE) %>%
    select(-Group)

####################
## Load genotypes ##
####################

rds.file <- sprintf("%s/augmented_data_big/ukb_gen_%s.rds", scratch, resolution)
if(file.exists(rds.file)){
    cat(sprintf("Found FBM in %s.\n", rds.file))
} else {
    cat(sprintf("Could not find FBM in %s.\n", rds.file))
    quit()
}

# Attach the "bigSNP" object in R session
ptm <- proc.time() # Start the clock!
cat("Attaching bigSNP object... ")
obj.bigSNP <- snp_attach(rds.file)
cat("done.\n")
proc.time() - ptm # Stop the clock

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
SNP <- obj.bigSNP$map$marker.ID

#####################
## Load phenotypes ##
#####################
cat("Reading phenotype file... ")

# Load list of subjects
Subjects <- obj.bigSNP$fam %>% as_tibble() %>% transmute(FID=family.ID, IID=sample.ID)

# Load response values
pheno.file <- sprintf("%s/simulations/phenotypes/%s_phenotypes.tab", scratch, experiment.name)
Phenotypes <- read_delim(pheno.file, delim="\t", col_types=cols())
keep.rows <- which(Phenotypes$FID %in% Subjects$FID)
Phenotypes <- Phenotypes %>% right_join(Subjects, by=c("FID", "IID"))

cat("done.\n")

##############################
## Load association results ##
##############################

# Read clumping file
clumped.file <- sprintf("%s/simulations/bolt_clumped/fold_%s/stats_%s_%s_%s.clumped",
                      scratch, fold, experiment.name, amplitude, clumping)
Clumped.raw <- read.table(clumped.file, header=TRUE) %>% as_tibble() %>%
    select(CHR, SNP, BP, P, SP2)

# Parse list of clumped discoveries
if(nrow(Clumped.raw)==0) {
    Clumped <- tibble()
} else {
    Variants.prefix <- Variants.original %>%
        mutate(SNP.sec = gsub(".A", "", SNP), SNP.sec = gsub(".B", "", SNP.sec), BP.sec=BP) %>%
        select(CHR, SNP, SNP.sec, BP.sec)
    # Extract lead SNPs
    Clumped.raw.lead <- Clumped.raw %>% select(CHR, BP, P) %>%
        inner_join(select(Variants.original, CHR, SNP, BP), by = c("CHR", "BP")) %>%
        mutate(SNP.lead=SNP, BP.lead=BP, SNP=SNP) %>%
        select(CHR, SNP.lead, BP.lead, P, SNP, BP)
    # Extract secondary SNPs
    Clumped.raw.secondary <- Clumped.raw %>%
        mutate(SNP.lead=SNP, SNP.sec = gsub("\\(1\\)", "", SP2)) %>%
        separate_rows(SNP.sec, sep=",") %>%
        select(CHR, BP, P, SNP.sec) %>%
        filter(SNP.sec!="NONE") %>%
        inner_join(select(Variants.original, CHR, SNP, BP), by = c("CHR", "BP")) %>%
        mutate(SNP.lead=SNP, BP.lead=BP) %>%
        select(CHR, SNP.lead, BP.lead, P, SNP.sec) %>%
        inner_join(Variants.prefix, by = c("CHR", "SNP.sec")) %>%
        mutate(BP=BP.sec) %>%
        select(CHR, SNP.lead, BP.lead, P, SNP, BP, P)
    # Combine lists of lead and secondary SNPs
    Clumped <- rbind(Clumped.raw.lead, Clumped.raw.secondary) %>% as_tibble() %>%
        arrange(CHR, BP.lead, BP)
    # Cross-reference with complete list of variants
    Clumped <- Clumped %>%
        inner_join(select(Variants.original,CHR,SNP,BP,Causal), by = c("CHR", "SNP", "BP"))
}

cat(sprintf("Loaded %d clumps containing %d causal variants\n",
            length(unique(Clumped$SNP.lead)), sum(Clumped$Causal)))

loci.intersect <- function(bp.min.x, bp.max.x, bp.min.y, bp.max.y, gap=1e5) {
    intersect <- as.logical((bp.min.y<=bp.max.x+gap) * (bp.max.y>=bp.min.x-gap))
    return(intersect)
}

if(nrow(Clumped)>0) {

    # Consolidate nearby clumps
    cat(sprintf("Consolidating nearby clumps... "))
     Selections.consolidated <- Clumped %>%
         group_by(CHR, SNP.lead, BP.lead) %>%
         summarise(P=min(P), BP.min=min(BP), BP.max=max(BP), Size=n(), Causal=any(Causal)) %>%
         mutate(BP.width=BP.max-BP.min) %>%
         ungroup() %>%
         consolidate_clumps(gap=1e5)

     # Cross-reference consolidated selections with clumps
     Clumped.consolidated <- full_join(Clumped %>% select(CHR, SNP, BP, P, Causal),
               Selections.consolidated %>% select(CHR, SNP.lead, BP.lead, BP.min, BP.max),
               by = c("CHR")
               ) %>%
         filter(BP>=BP.min, BP<=BP.max)

    cat("done.\n")
    cat(sprintf("There are now %d clumps containing %.1f variants each, on average.\n",
                length(unique(Clumped.consolidated$SNP.lead)),
                Clumped.consolidated %>% group_by(SNP.lead) %>% summarise(N=n()) %>%
                summarise(N=mean(N)) %>% as.numeric()))

     # Sanity check
     stopifnot(nrow(Clumped.consolidated)==nrow(Clumped))
     stopifnot(length(unique(Clumped.consolidated$SNP.lead))==length(unique(Selections.consolidated$SNP.lead)))

} else {
    Clumped.consolidated <- Clumped
}

#####################################
## Perform fine mapping with SUSIE ##
#####################################

coverage <- 0.9
min_abs_corr <- 0.1
correct.bias <- TRUE

if(nrow(Clumped.consolidated)>0) {
    Clumped.consolidated <- Clumped.consolidated %>% select(CHR, BP, SNP.lead, BP.lead) %>%
        inner_join(Variants.original %>% select(CHR, SNP, BP), by = c("CHR", "BP")) %>%
        select(CHR, SNP.lead, BP.lead, SNP, BP, everything())

    if(correct.bias) {
        # Insert all the missing variants in each clump, to correct for selection bias
        Clumped.full <- lapply(unique(Clumped.consolidated$SNP.lead), function(clump.snp.lead) {
            Clump <- Clumped.consolidated %>% filter(SNP.lead==clump.snp.lead)
            # Add the missing variants to the clump, to mitigate selection bias
            clump.chr <- Clump$CHR[1]
            clump.bp.min <- min(Clump$BP)
            clump.bp.max <- max(Clump$BP)
            clump.snp.lead <- Clump$SNP.lead[1]
            clump.bp.lead <- Clump$BP.lead[1]
            Clump <- Variants.original %>%
                filter(CHR==clump.chr, BP>=clump.bp.min, BP<=clump.bp.max) %>%
                mutate(SNP.lead=clump.snp.lead, BP.lead=clump.bp.lead) %>%
                select(CHR, SNP.lead, BP.lead, SNP, BP)
            return(Clump)
        })
        Clumped.full <- do.call("rbind", Clumped.full)
        Clumped.consolidated <- Clumped.full
    }

    cat(sprintf("Running SUSIE on %d clumps with coverage %.2f and min_abs_corr %.2f\n",
                length(unique(Clumped.consolidated$SNP.lead)), coverage, min_abs_corr))
    cat(sprintf("The %d clumps containing %.1f variants each, on average.\n",
                length(unique(Clumped.consolidated$SNP.lead)),
                Clumped.consolidated %>% group_by(SNP.lead) %>% summarise(N=n()) %>%
                summarise(N=mean(N)) %>% as.numeric()))

    # Apply SUSIE clump by clump
    n.clumps <- length(unique(Clumped.consolidated$SNP.lead))
    Clumped.fine <- lapply(1:n.clumps, function(clump.idx) {
        clump.snp.lead <- unique(Clumped.consolidated$SNP.lead)[clump.idx]
        Clump <- Clumped.consolidated %>% filter(SNP.lead==clump.snp.lead)
        cat(sprintf("Running SUSIE on clump %d of %d (size %d) ...\n", clump.idx, n.clumps, nrow(Clump)))
        # Apply fine mapping to this clump
        Clump.fine <- fine.mapping(Clump, obj.bigSNP, Phenotypes, pheno.name,
                                   coverage=coverage, min_abs_corr=min_abs_corr, simplify=TRUE)
        # Add ground truth
        Clump.fine <- Clump.fine %>%
            inner_join(select(Variants.original, CHR, SNP, Causal), by = c("CHR", "SNP"))
        # Count false discoveries on this clump
        n.discoveries <- length(unique(Clump.fine$SNP.lead))
        n.discoveries.false <- Clump.fine %>% group_by(SNP.lead) %>% summarise(Causal=any(Causal)) %>%
                    summarise(False=sum(!Causal)) %>% as.numeric()
        if(n.discoveries.false>0) {
            cat(sprintf("False discoveries made by SUSIE: %d out of %d \n",
                        n.discoveries.false, n.discoveries))
        }
        return(Clump.fine)
    })
    Clumped.fine <- do.call("rbind", Clumped.fine)
        
    # Print summary of results
    cat("Summary of discoveries made by SUSIE:\n")
    Clumped.fine %>%
        group_by(CHR, SNP.lead) %>% summarise(Causal=any(Causal), Size=n()) %>%
        ungroup() %>% summarise(Size=mean(Size), True=sum(Causal), FDP=mean(!Causal)) %>%
        print()
    
} else {
    Clumped.fine <- tibble()
}

# Save results
if(correct.bias) {
    out.file <- sprintf("%s/simulations/susie/fold_%s/%s_%s_%s_%s_nobias.txt",
                        scratch, fold, experiment.name, amplitude, clumping, coverage)
} else {
    out.file <- sprintf("%s/simulations/susie/fold_%s/%s_%s_%s_%s.txt",
                        scratch, fold, experiment.name, amplitude, clumping, coverage)
}
Clumped.fine %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written on:\n %s\n", out.file))

# Test
if(FALSE) {

    # With bias correction
    clump.snp.lead <- unique(Clumped.full$SNP.lead)[1]
    Clump.in <- Clumped.full %>% filter(SNP.lead==clump.snp.lead)
    Clump.fine <- Clump.in %>%
        fine.mapping(obj.bigSNP, Phenotypes, pheno.name,
                     coverage=coverage, min_abs_corr=min_abs_corr, simplify=TRUE)
    Clump.fine %>% inner_join(select(Variants.original, CHR, SNP, Causal), by = c("CHR", "SNP"))

    # Without bias correction
    clump.snp.lead <- unique(Clumped.consolidated$SNP.lead)[1]
    Clump.in <- Clumped.consolidated %>% filter(SNP.lead==clump.snp.lead)
    Clump.fine <- Clump.in %>%
        fine.mapping(obj.bigSNP, Phenotypes, pheno.name,
                     coverage=coverage, min_abs_corr=min_abs_corr, simplify=TRUE)
    Clump.fine %>% inner_join(select(Variants.original, CHR, SNP, Causal), by = c("CHR", "SNP"))

}
