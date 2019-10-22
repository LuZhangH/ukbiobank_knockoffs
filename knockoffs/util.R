load.unphased <- function(chr, select.subjects=NULL, select.snps=NULL) {
    # Names of files where the original data is stored
    data.dir <- paste("/scratch/PI/candes/ukbiobank/genotypes/",sep="")
    bed.file <- paste(data.dir, "ukb_gen_chr", chr, ".bed", sep = "")
    bim.file<- paste(data.dir, "ukb_gen_chr", chr, ".bim", sep = "")
    fam.file <- paste(data.dir, "ukb_gen_chr", chr, ".fam", sep = "")

    # Select relevant SNPs
    if(!is.null(select.snps)) {
        Variants <- read_delim(bim.file, delim="\t", col_names=c("Chr","Variant","X1","Position","A1","A2"),
                               progress = FALSE, col_types = cols())
        select.snps <- which(Variants$Variant %in% select.snps)
    }
    
    # Load original genotypes
    Data <- read.plink(bed.file, bim.file, fam.file, select.snps=select.snps)
    Data$fam <- as.tibble(Data$fam)
    Data$map <- as.tibble(Data$map)
    colnames(Data$map) <- c("Chr","Variant","X1","Position","A1","A2")

    # Select relevant subjects
    if(!is.null(select.subjects)) {
        select.subjects <- which(Data$fam$pedigree %in% select.subjects)
        Data$fam <- Data$fam[select.subjects,]
        Data$genotypes <- Data$genotypes[select.subjects,]
    }
    Data$genotypes <- as(Data$genotypes, "numeric")
    
    return(Data)
}

load.phased <- function(fp.dir, chr, n.skip=0, n.load=-1, select.snps=NULL) {
    Data = c()

    # Load haplotype data
    inp.file <- paste(fp.dir, "input/ukb_hap_chr", chr, ".inp", sep = "")
    if( n.load<0 ) {
        n.load <- as.numeric(read_lines(inp.file, n_max=1, skip=0))
    }
    p <- as.numeric(read_lines(inp.file, n_max=1, skip=1))
    Data$haplotypes <- read_fwf(inp.file, skip=3+3*n.skip, n_max=2*n.load, fwf_widths(rep(1,p)),
                               comment="#", col_types=cols(.default = col_integer()), progress = FALSE)
    Data$haplotypes <- as.matrix(Data$haplotypes)

    # Load subject metadata
    fam.file <- paste(fp.dir, "input/ukb_hap_chr",chr,".fam", sep = "")
    Data$fam <- read_delim(fam.file, delim=" ", skip=n.skip, n_max=n.load, progress = FALSE,
                           col_types=cols(.default = col_integer()),
                           col_names=c("pedigree", "member","X1","X2","X3","X4"))
    # Load SNP metadata
    positions <- scan(inp.file, what=integer(), skip=2, nlines=1, sep=" ", na.strings=c("P"), quiet=TRUE)
    Data$map  <- tibble("Position"=positions[2:(length(positions)-1)])

    # Select relevant SNPs
    if(!is.null(select.snps)) {
        select.snps <- which(Data$map$Position %in% select.snps)
        Data$map <- Data$map[select.snps,]
        Data$haplotypes <- Data$haplotypes[,select.snps]
    }

    return(Data)
}

determine.encoding <- function(fp.dir, chr, n.skip, n.load) {
    # Load block of phased data
    dat.phased <- load.phased(fp.dir, chr, n.skip=n.skip, n.load=n.load,
                              select.snps=Variants$Position)
    dat.phased$genotypes <- dat.phased$haplotypes[seq(1,nrow(dat.phased$haplotypes),by=2),] +
        dat.phased$haplotypes[seq(2,nrow(dat.phased$haplotypes),by=2),]
    # Load block of unphased data
    dat.unphased <- load.unphased(chr, select.snps=Variants$Variant, select.subjects=dat.phased$fam$pedigree)
    # Match values
    flip.sign <- which(sapply(1:ncol(dat.unphased$genotypes), function(j) {
        not.missing <- which(!is.na(dat.unphased$genotypes[,j]))
        sum(dat.unphased$genotypes[not.missing,j] != dat.phased$genotypes[not.missing,j]) != 0
    }))
    return(flip.sign)
}

