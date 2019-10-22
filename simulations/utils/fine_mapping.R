#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

source("../utils/util.R")
suppressMessages(library(tidyverse))
suppressMessages(library(susieR))
suppressMessages(require(igraph))
suppressMessages(library(doParallel))

simplify.nested.cs <- function(sets.cs) {
    # Remove nested sets
    if(is.null(sets.cs)) {
        return(NULL)
    }
    # Sort each set
    sets.cs <- lapply(sets.cs, function(set) sort(set))
    # Remove identical sets
    sets.cs <- unique(sets.cs)
    # Count number of distinct sets
    n.sets <- length(sets.cs)

    # Check whether any sets are contained in other sets
    contained <- sapply(1:n.sets, function(i) sapply(1:n.sets, function(j) {
        if(i==j) {
            return(0)
        }
        if(length(setdiff(sets.cs[[j]],sets.cs[[i]]))==0) {
            return(1)
        }
        return(0)
    }))
    contained <- matrix(contained, n.sets, n.sets)

    # Remove sets that contain another set
    sets.simplified <- unique(lapply(seq(n.sets), function(i) {
        remove.idx <- which(contained[,i]==1)
        if(sum(remove.idx)>0){
            return(NULL)
        } else {
            return(sort(sets.cs[[i]]))
        }

    }))
    sets.is.null <- sapply(sets.simplified, function(set) is.null(set))
    sets.simplified <- sets.simplified[!sets.is.null]
    return(sets.simplified)
}

merge.overlapping.cs <- function(sets.cs) {
    # Merge sets with partial intersection
    if(is.null(sets.cs)) {
        return(NULL)
    }
    # Check pairwise intersections between all reported sets
    n.sets <- length(sets.cs)
    intersections <- sapply(1:n.sets, function(i) sapply(1:n.sets, function(j) {
        if(i==j) {
            return(1)
        }
        if(length(intersect(sets.cs[[i]],sets.cs[[j]]))>0) {
            return(1)
        }
        return(0)
    }))
    intersections <- matrix(intersections, n.sets, n.sets)

    # Extract connected components
    g <- igraph::graph.adjacency(intersections)
    g.clu <- igraph::components(g)
    g.groups <- igraph::groups(g.clu)
    sets.distinct <- lapply(g.groups, function(group) {
        sort(unique(unlist(lapply(group, function(g) sets.cs[[g]]))))
    })
    names(sets.distinct) <- NULL
    return(sets.distinct)
}

comma_sep = function(x) {
    x = strsplit(as.character(x), ",")
    unlist(lapply(x, paste, collapse = ','))
}

fine.mapping <- function(Clump.inp, obj.bigSNP, Phenotypes, pheno.name,
                         coverage=0.9, min_abs_corr=0.1, simplify=FALSE) {
  
    # Make sure that the rows of the genotypes match the rows of the phenotypes
    Phenotypes.full <- Phenotypes %>%
        right_join(transmute(obj.bigSNP$fam, FID=family.ID, IID=sample.ID), by = c("FID", "IID"))

    # Compute row and column indices
    ind.row <- which(obj.bigSNP$fam$family.ID %in% Phenotypes.full$FID)
    ind.col <- which(obj.bigSNP$map$marker.ID %in% Clump.inp$SNP)

    # Extract relevant genotypes and phenotypes
    X <- obj.bigSNP$genotypes[ind.row, ind.col, drop=FALSE]
    y <- Phenotypes.full[[pheno.name]][ind.row]

    susie.fitted <- susieR::susie(X, y, L=10, estimate_residual_variance=TRUE, scaled_prior_variance=0.1)
    
    # Extract credible sets
    susie.sets <- susieR::susie_get_cs(susie.fitted, X=X, coverage=coverage, min_abs_corr=min_abs_corr)
    sets.cs <- susie.sets$cs
    cat(sprintf("SUSIE returned %d credible sets for %d variables.\n", length(sets.cs), nrow(Clump.inp)))
    names(sets.cs) <- NULL

    # Simplify overlapping sets
    if(simplify==TRUE) {
        sets.cs <- simplify.nested.cs(sets.cs)
        sets.cs <- merge.overlapping.cs(sets.cs)
        cat(sprintf("The credible sets have been simplified: %d remaining.\n", length(sets.cs)))
    } else {
        cat(sprintf("The credible sets have not been simplified.\n"))
    }

    # Summarise results
    if(is.null(sets.cs)) {
        Clump.fine <- Clump.inp %>% select(CHR, SNP.lead, BP.lead, SNP) %>% head(0)
    } else {
        sets.leads <- sapply(sets.cs, function(variants) {variants[1]})
        sets.snps <- lapply(sets.cs, function(variants) {Clump.inp$SNP[variants]})

        Susie.sets <- tibble(CHR=Clump.inp$CHR[1], SNP.lead=Clump.inp$SNP[sets.leads],
                             BP.lead=Clump.inp$BP[sets.leads], SNPs=sets.snps)

        Clump.fine <- Susie.sets %>% mutate(SNPs = comma_sep(SNPs)) %>%
            separate_rows(SNPs, sep=",") %>% mutate(SNP = SNPs) %>% select(-SNPs) %>%
            mutate(SNP=gsub(" ", "", SNP, fixed = TRUE)) %>%
            mutate(SNP=gsub("\"", "", SNP), SNP=gsub("c\\(", "", SNP), SNP=gsub(")", "", SNP))
    }
    return(Clump.fine)
}

fine.mapping.caviar <- function(Clumped, Association, LD.full) {
    snp.lead.list <- unique(Clumped$SNP.lead)

    tmp.dir <- tempdir()

    #Posterior <- lapply(1:length(snp.lead.list), function(clump.num) {
    Posterior <- foreach(clump.num = 1:length(snp.lead.list), .combine = 'rbind') %dopar% {
        cat(sprintf("Performing fine mapping on clump %d of %d ...\n", clump.num, length(snp.lead.list)))

        snp.lead <- snp.lead.list[clump.num]

        Clump <- Clumped %>% filter(SNP.lead==snp.lead) %>% select(CHR, SNP.lead, BP.lead, SNP, BP) %>%
            left_join(select(Association, CHR, BP, STAT, P), by=c("CHR", "BP")) %>%
            arrange(BP)

        # Convert LD to a symmetric square matrix
        LD <- filter(LD.full, BP_A %in% Clump$BP, BP_B %in% Clump$BP)
        LD.variants <- unique(c(LD$SNP_A, LD$SNP_B))
        LD.positions <- unique(c(LD$BP_A, LD$BP_B))

        Sigma <- ld.to.mat(LD, "R", Clump$BP)

        ##############################
        ## Prepare input for CAVIAR ##
        ##############################

        caviar.inp.Z <- sprintf("%s/Z_%s.txt", tmp.dir, clump.num)
        Clump %>% select(SNP, STAT) %>% write_delim(caviar.inp.Z, delim="\t", col_names=FALSE)

        caviar.inp.LD <- sprintf("%s/LD_%s.txt", tmp.dir, clump.num)
        Sigma <- Sigma %>% as.matrix()
        dimnames(Sigma)=NULL
        Sigma %>% write.matrix(caviar.inp.LD, sep="\t")

        caviar.out <- sprintf("%s/caviar_out_%s", tmp.dir, clump.num)
        caviar.command <- sprintf("ml gsl; CAVIAR -c 2 -r 0.9 -l %s -z %s -o %s",
                                  caviar.inp.LD, caviar.inp.Z, caviar.out)
        system(caviar.command)

        # Load CAVIAR output
        Posterior.clump <- read.table(sprintf("%s_post", caviar.out), header=T,
                                col.names=c("SNP", "P.set", "P.causal")) %>%
            as_tibble() %>% mutate(SNP=as.character(SNP)) %>%
            left_join(Clump, by=c("SNP")) %>%
            select(CHR, SNP.lead, BP.lead, SNP, BP, STAT, P, P.set, P.causal)

        # Remove CAVIAR files
        clean.command <- sprintf("rm %s %s %s_post %s.log %s_set",
                                 caviar.inp.LD, caviar.inp.Z, caviar.out, caviar.out, caviar.out)
        system(clean.command)

        return(Posterior.clump)
    }
    #Posterior <- do.call("rbind", Posterior)
    return(Posterior)
}
