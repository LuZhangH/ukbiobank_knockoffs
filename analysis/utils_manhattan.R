knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset>1 | offset<0) {
    stop('Input offset must be between 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
                 
pre_process <- function(gwas) {
    gwas.don <- gwas %>% 
        # Compute chromosome size
        group_by(CHR) %>% 
        summarise(chr_len=max(BP)) %>% 

        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%

        # Add this info to the initial dataset
        left_join(gwas, ., by=c("CHR"="CHR")) %>%

        # Add a cumulative position of each SNP
        arrange(CHR, BP) %>%
        mutate( BPcum=BP+tot)

    return(gwas.don)
}

manhattan_simple <- function(gwas.don, axisdf, limit.left, limit.right, ytrans="log10") {
    # Make LMM Manhattan plot
    y.max <- 1.1 * max(-log10(pmax(1e-300,gwas.don$P)))
    p.manhattan <- gwas.don %>%
        mutate(P=pmax(1e-300,P)) %>%
        filter(P<1e-3) %>%
        ggplot(aes(x=BPcum, y=-log10(P))) +

        # Show significance threshold
        geom_hline(yintercept=-log10(5e-8), linetype="dashed", color = "red") +
    
        # Show all points
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
        scale_color_manual(values = rep(c("darkgrey", "black"), 22 )) +

        # Custom axes:
        scale_x_continuous(label = axisdf$CHR, breaks=axisdf$center, limits=c(limit.left,limit.right)) +
        scale_y_continuous(trans=ytrans, expand = c(0,0), limits=c(3,y.max)) +     # remove space between plot area and x axis
        xlab("Chromosome") + ylab("Importance") +
    
        # Custom the theme:
        theme_bw() +
        theme(legend.position="none", panel.border = element_blank(), 
              panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
             )

    return(p.manhattan)
}

plot_manhattan <- function(gwas, ytrans="log10") {
    
    don <- pre_process(gwas)
    
    # Compute plot limits
    limit.left <- min(don$BPcum)
    limit.right <- max(don$BPcum)
    
    # Find centers of each chromosome
    axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    # Make plot
    p.manhattan <- manhattan_simple(don, axisdf, limit.left, limit.right, ytrans=ytrans)
    
    return(p.manhattan)
}

plot_manhattan_knockoffs <- function(LMM, Knockoffs, ytrans="log10") {
    
    # Create chromosome blocks
    LMM.don <- pre_process(LMM)
    
    # Compute plot limits
    limit.left <- min(LMM.don$BPcum)
    limit.right <- max(LMM.don$BPcum)
    
    # Find centers of each chromosome
    axisdf = LMM.don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    # Make LMM Manhattan plot
    p.manhattan <- manhattan_simple(LMM.don, axisdf, limit.left, limit.right, ytrans=ytrans)
    
    # Transform knockoffs results into pseudo p-values
    W.thresh <- knockoff.threshold(Knockoffs$W)
    if(W.thresh==Inf) {
        W.thresh <- max(Stats$W)*1.1
    }
    W.scale <- -log(5e-8)/W.thresh
    Knockoffs <- Knockoffs %>% 
        filter(W>0) %>% 
        mutate(SNP=SNP.lead, BP=BP.lead, P.W=exp(-W*W.scale)) %>%
        left_join(LMM, ., by = c("SNP", "CHR", "BP")) %>%
        mutate(P=P.W)
    
    p.is.na <- which(is.na(Knockoffs$P))
    Knockoffs$P[p.is.na] <- runif(length(p.is.na))
    
    # Create chromosome blocks
    Knockoffs.don <- pre_process(Knockoffs)
    
    # Find centers of each chromosome
    axisdf = Knockoffs.don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    # Make Knockoffs Manhattan plot
    p.knockoffs <- manhattan_simple(Knockoffs.don, axisdf, limit.left, limit.right, ytrans=ytrans)
    
    plot_grid(p.manhattan, p.knockoffs, ncol=1, align="v", axis="tblr") 
}