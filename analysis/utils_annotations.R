place_gene_segments <- function(Genes.window, verbose=FALSE) {
    
    # Define function that checks whether a new gene would fit in an existing row
    gene_fits <- function(txStart, txEnd, row, gap=5e4) {
        if(length(row)==0) {
            return(TRUE)
        }
        for(segment in row) {
            if(txStart<=segment[2]+gap && txEnd>=segment[1]-gap){
                return(FALSE)
            }
        }
        return(TRUE)
    }
    
    # Sort genes by length
    Genes.window <- Genes.window %>% mutate(txWidth = txEnd-txStart)
    gene.order <- order(Genes.window$txWidth,decreasing=T)
    Genes.window <- Genes.window %>% arrange(desc(txWidth))
    
    # Count number of genes
    n.genes <- nrow(Genes.window)
    gene.heights <- rep(NA, n.genes)
    gene.rows <- list()
    
    n.rows <- 0
    
    for(j in 1:n.genes) {
        txStart <- Genes.window$txStart[j]
        txEnd <- Genes.window$txEnd[j]

        if(length(gene.rows)==0) {
            # Add the first gene to the first row
            gene.rows[[1]] <- list()
            gene.rows[[1]][[1]] <- c(txStart,txEnd)
            gene.heights[j] <- 1
            n.rows <- n.rows + 1
        } else {
            would.fit <- sapply(1:length(gene.rows), function(g) gene_fits(txStart, txEnd, gene.rows[[g]]))
            if(any(would.fit)) {
                # Add gene to existing row
                k <- min(which(would.fit))
                #print(length(gene.rows[[k]]))
                row.length <- length(gene.rows[[k]])
                if(verbose) cat(sprintf("Appending (%d,%d) to row %d \n",txStart,txEnd,k))
                gene.rows[[k]][[row.length+1]] <- c(txStart,txEnd)
                gene.heights[j] <- k
            } else {
                # Add gene to new row
                n.rows <- n.rows + 1
                if(verbose) cat(sprintf("No room to append (%d,%d). Creating new row %d. \n",txStart,txEnd,n.rows))
                gene.rows[[n.rows]] <- list()
                gene.rows[[n.rows]][[1]] <- c(txStart,txEnd)
                gene.heights[j] <- n.rows
            }
        }
    }

    # Assign heights to each gene
    Genes.window$Height <- -(gene.heights-0.5)
                                
    # Reorganise rows in the original order
    Genes.window[gene.order,] <- Genes.window
    
    # Return modified data
    return(Genes.window)
}
                                
plot_sears_tower <- function(window.chr, window.left, window.right, 
                             Discoveries, LMM, Clumped, Annotations.func, Exons.canonical,
                             plot.file=NULL) {

   p.significant <- 5e-8
    
   # Extract color map
    annotation.color.map <- Annotations.func %>% group_by(name, itemColor) %>% summarise() %>% 
        ungroup() %>%
        mutate(name.num=parse_number(as.character(name))) %>% 
        mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
        mutate(label=as.factor(label)) %>%
        arrange(name.num)
    
    # Convert names to factors according to color maps
    Annotations.func <- Annotations.func %>%
        mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
        mutate(label=factor(label, levels=annotation.color.map$label, labels=annotation.color.map$label))
    
    # Extract knockoff discoveries within this window
    Knockoffs.window <- Discoveries %>% filter(Method=="Knockoffs") %>%
        filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
    cat(sprintf("There are %d knockoff discoveries within this window.\n", nrow(Knockoffs.window)))

    # Update window limits
    window.left <- min(window.left, min(Knockoffs.window$BP.min))
    window.right <- max(window.right, max(Knockoffs.window$BP.max))

    # Extract LMM pvalues within this window
    LMM.clumped.window <- Clumped %>% filter(CHR==window.chr, BP<=window.right, BP>=window.left)
    LMM.window <- LMM %>% filter(CHR==window.chr, BP<=window.right, BP>=window.left) %>%
        left_join(LMM.clumped.window, by = c("SNP", "CHR", "BP")) %>% 
        mutate(BP.lead=factor(BP.lead))
    
    cat(sprintf("There are %d LMM pvalues within this window, %d of which are significant.\n", 
                nrow(LMM.window), sum(LMM.window$P<p.significant)))

    # Select exons within this windows
    Exons.window <- Exons.canonical %>% 
        filter(chrom==window.chr, txStart<=window.right, txEnd>=window.left)
    cat(sprintf("There are %d exons within this window, divided into %d genes.\n", 
                nrow(Exons.window), length(unique(Exons.window$name2))))

    # Select functional annotations within this window
    Functional.window <- Annotations.func %>%
        filter(chrom==window.chr, chromStart<=window.right, chromEnd>=window.left)
    cat(sprintf("There are %d functional annotations within this window.\n", 
                nrow(Functional.window)))

    # Significance level for pvalues
    Window.nominal <- LMM.window %>% mutate(Importance = -log10(p.significant))

    # Minimal theme
    theme_minimal <- theme_bw() + 
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.border=element_blank()
             )
    
    # Manhattan plot
    LMM.window$BP.lead <- round((LMM.window$BP.lead %>% as.character %>% parse_number)/1e6,3) %>% factor
    
    if(all(is.na(LMM.window$BP.lead))) {
        p.manhattan <- LMM.window %>%
            mutate(P=pmax(1e-300,P)) %>%
            ggplot(aes(x=BP, y=-log10(P))) +
            geom_point(color="black", alpha=0.25)
    } else {
        p.manhattan <- LMM.window %>%
            mutate(P=pmax(1e-300,P)) %>%
            ggplot(aes(x=BP, y=-log10(P), color=BP.lead, alpha=is.na(BP.lead))) +
            geom_point() +
            scale_colour_discrete(na.value = "black", name="Clump lead", guide=FALSE)
    }
    p.manhattan <- p.manhattan +
        geom_hline(data=Window.nominal, aes(yintercept=Importance), linetype="dashed", color = "black") +
        scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 1), guide=FALSE) +
        scale_x_continuous(limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
        scale_y_continuous(limits=c(0,300), trans="sqrt") +
        theme_bw() +
        xlab(sprintf("Chromosome %d (Mb)", window.chr)) + ylab(TeX("-log_{10}(p)")) +
        theme(panel.border = element_blank(), 
              plot.title = element_text(size = 9,face="bold"), 
              axis.title = element_text(size=8)) +
        ggtitle(sprintf("Manhattan plot (%s)", phenotype))

    # Plot clumped LMM discoveries
    LMM.clumps.window <- LMM.clumped.window %>% group_by(CHR, SNP.lead, BP.lead) %>%
        summarise(BP.min=min(BP), BP.max=max(BP)) %>%
        ungroup() %>%
        mutate(BP.lead=factor(BP.lead))

    clump.list <- levels(LMM.clumps.window$BP.lead)
    clump.heights <- seq(length(clump.list)) %>% rev
    names(clump.heights) <- clump.list
    
    p.clumped <- LMM.clumps.window %>%
        mutate(Height=clump.heights[BP.lead]) %>%
        ggplot() +
        geom_segment(aes(x=pmax(BP.min,window.left), y=Height, xend=pmin(BP.max,window.right), yend=Height, color=BP.lead)) +
        geom_segment(aes(x=BP.min, y=Height-0.4, xend=BP.min, yend=Height+0.4, color=BP.lead)) +
        geom_segment(aes(x=BP.max, y=Height-0.4, xend=BP.max, yend=Height+0.4, color=BP.lead)) +
        scale_colour_discrete(na.value = "black", guide=FALSE) +
        scale_x_continuous(limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
        scale_y_continuous(breaks = NULL) +
        ylab("") + xlab("") +
        theme_void() +
        ggtitle("Clumped discoveries") + 
        theme(plot.title = element_text(size = 9,face="bold"))
    
    # Plot knockoff discoveries
    resolution.list <- levels(Knockoffs.window$Resolution)
    resolution.heights <- seq(length(resolution.list)) %>% rev
    names(resolution.heights) <- resolution.list
    
    p.knockoffs <- Knockoffs.window %>% 
        mutate(Height=resolution.heights[Resolution]) %>%
        ggplot() +
        geom_rect(aes(xmin=pmax(BP.min,window.left), xmax=pmin(BP.max,window.right), ymin=Height-0.5, ymax=Height+0.5), 
                  alpha=0.5, fill="black", color="black") +
        ylab("") + xlab("") +
        scale_x_continuous(limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
        scale_y_continuous(limits=c(0.5,max(resolution.heights)-0.5), breaks = NULL) +
        theme_minimal +
        ggtitle("Discoveries with knockoffs") + 
        theme(plot.title = element_text(size = 9,face="bold"))
    
    # Plot functional annotations
    myColors <- annotation.color.map$itemColor
    names(myColors) <- annotation.color.map$label
    
    p.functional <- Functional.window %>%
        mutate(chromStart=pmax(chromStart, window.left), chromEnd=pmin(chromEnd, window.right)) %>%
        ggplot() +
        geom_rect(aes(xmin=chromStart, xmax=chromEnd, ymin=0.5, ymax=1.5, fill=label)) +
        ylab("") + xlab("") +
        scale_x_continuous(limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
        scale_color_manual(values=myColors, guide=FALSE) +
        scale_fill_manual(values=myColors, name="Function") +
        theme_void() +
        theme(legend.key.size = unit(0.5, "cm"), legend.title=element_text(size=8,face="bold"), 
              legend.text=element_text(size=8),
              plot.title = element_text(size = 9,face="bold")) +
        ggtitle("Functional annotations")
    
    Genes.window <- Exons.window %>% group_by(name, name2, strand) %>%
        summarise(txStart=min(txStart), txEnd=max(txEnd)) %>%
        mutate(txStart=max(txStart, window.left), txEnd=min(txEnd, window.right)) %>%
        place_gene_segments() %>%
        mutate(txCenter=(txStart+txEnd)/2, Height=Height/2)
    
    # Plot exons and genes
    p.genes <- Exons.window %>%
        filter(exonStarts>=window.left, exonEnds<=window.right) %>%
        inner_join(Genes.window %>% select(name, name2, Height), by = c("name", "name2")) %>%
        ggplot() +
        geom_rect(aes(xmin=exonStarts, xmax=exonEnds, ymin=Height-0.1, ymax=Height+0.1), 
                  alpha=1, color="black", fill="black") +
        geom_segment(data=Genes.window, aes(x=txStart, y=Height, xend=txEnd, yend=Height, group=name2), color="black") +
        geom_label_repel(data=Genes.window, aes(x=txCenter, y=Height, label=paste(name2,strand,sep=" ")),
                         size=1.5, direction="y", box.padding=0.35, point.padding=0.35, label.padding=0.1, force=2, 
                         segment.color = 'grey50', segment.alpha=0.5, seed=2019) +
        ylab("") + xlab("") +
        scale_x_continuous(limits=c(window.left,window.right), labels=function(x) {x*1e-6}) +
        scale_y_continuous(breaks = NULL) +
        theme_void() +
        ggtitle("Gene annotations") + 
        theme(plot.title = element_text(size = 9,face="bold"))
    
    #rel.height.clumps <- 0.1 + 0.1*length(clump.heights)
    #rel.height.knockoffs <- 0.25 + 0.1*length(unique(Knockoffs.window$Resolution))
    #rel.height.functional <- 0.1 + 0.5*length(unique(Genes.window$Height))
       
    #p <- plot_grid(p.manhattan, p.clumped, p.knockoffs, p.functional, p.genes, ncol=1, 
    #               align="v", axis="tblr", 
    #               rel_heights=c(2,rel.height.clumps,rel.height.knockoffs,0.4,rel.height.functional))
    
    height.manhattan <- 1.25
    height.clumps <- 0.2
    height.knockoffs <- 1.25
    height.functional <- 0.4
    #height.genes <- 0.35 + 0.35*length(unique(Genes.window$Height))
    height.genes <- 2
    
    heights <- c(height.manhattan, height.clumps, height.knockoffs, height.functional, height.genes)
    
    p.towers <- egg::ggarrange(p.manhattan, p.clumped, p.knockoffs, p.functional, p.genes, ncol=1, draw=F,
                       heights=heights)
   
    if(!is.null(plot.file)) {
        ggsave(plot.file, p.towers, width=7, height=sum(heights))
    }
    
    return(p.towers)
}