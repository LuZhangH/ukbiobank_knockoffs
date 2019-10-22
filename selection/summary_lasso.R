library(tidyverse)
library(gridExtra)
library(knockoff)

phenotype <- "bmi"
clumping <- "Radj50"
K <- 50
tmp.dir <- "/scratch/PI/candes/ukbiobank_tmp"
regression <- "linear"
offset <- 1

# Load feature importance measures
res.file <- sprintf("%s/association/%s/%s_K%d/lasso.txt", tmp.dir, phenotype, clumping, K)
Lasso <- read_delim(res.file, delim=" ", col_types=cols(), guess_max=1000000)

# Compute knockoff statistics
W.stats <- function(importance, knockoff) {
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w = z-zk                
}
Results.knockoffs <- Lasso %>% filter(!is.na(CHR)) %>%
    group_by(CHR, Group) %>%
    summarize(W = W.stats(Importance,Knockoff), SNP=SNP[1], BP=BP[1],
              Sign.W = factor(sign(W), levels=c(-1,0,1))) %>%
    ungroup() %>%
    arrange(desc(abs(W))) %>%
    select(CHR, SNP, BP, Group, W, Sign.W)

# Apply the knockoff filter
W.thres <- knockoff.threshold(Results.knockoffs$W, offset=0)
Selected.knockoffs <- Results.knockoffs %>% filter(W >= W.thres)

# Load marginal p-values
load.results <- function(chr) {
    cat(sprintf("Loading results for chromosome %d... ", chr))
    # Load p-values
    pvals.file <- sprintf("%s/association/%s/%s_K%d/ukb_chr%d.assoc.%s",
                          tmp.dir, phenotype, clumping, K, chr, regression)
    Results.raw <- read_table(pvals.file, guess_max=21474836, col_types=cols(), progress=FALSE)
    Results <- Results.raw %>% filter(TEST=="ADD") %>%
        separate(SNP, c("SNP", "Knockoff"), sep="\\.") %>%
        select(CHR, SNP, Knockoff, BP, P)

    # Load identity of knockoffs
    key.file <- sprintf("%s/knockoffs/%s_K%d/ukb_gen_chr%d.key", tmp.dir, clumping, K, chr)
    Legend <- read_delim(key.file, delim=" ", col_types=cols(), progress=FALSE)
    Results <- Results %>% mutate(Knockoff=Legend$Knockoff)

    # Load group information
    grp.file <- sprintf("%s/augmented_data/%s_K%d/ukb_gen_chr%d.grp", tmp.dir, clumping, K, chr)
    Groups <- read_delim(grp.file, delim=" ", col_types=cols(), progress=FALSE)
    Groups <- Groups %>% mutate(SNP=Variant, BP=Position, CHR=Chr) %>% select(CHR, SNP, BP, Group)
    Results <- inner_join(Results, Groups, by=c("CHR", "SNP", "BP"))

    cat("done.\n")
    return(Results)
}

Results.marginal <- tibble()
for(chr in seq(1,22)) {
    Results.chr <- load.results(chr)
    Results.marginal <- rbind(Results.marginal, Results.chr)
}

# Select groups with marginal testing
Selected.marginal <- Results.marginal %>% filter(Knockoff==FALSE) %>%
    group_by(CHR, Group) %>%
    summarise(P=min(P), SNP=SNP[1], BP=round(mean(BP)), Group.size=n()) %>%
    select(CHR, SNP, BP, P, Group, Group.size) %>%
    mutate(Importance = -log10(P)) %>%
    arrange(desc(Importance)) %>%
    filter(P < 5e-8)

# Compare findings
plot.manhattan <- function(chr) {
    Results.chr <- filter(Results.marginal,CHR==chr)
    x.min <- min(Results.chr$BP)
    x.max <- max(Results.chr$BP)

    # Combine discoveries
    Selected.chr <- full_join(Selected.knockoffs %>% filter(CHR==chr) %>%
                          mutate(Method = "Knockoffs", Threshold=W.thres, Importance=W) %>%
                          select(c("CHR", "SNP", "BP", "Group", "Importance", "Threshold", "Method")),
                          Selected.marginal %>% filter(CHR==chr) %>%
                          mutate(Method = "Marginal", Threshold=-log10(5e-8)) %>%
                                 select(c("CHR", "SNP", "BP", "Group", "Importance", "Threshold", "Method"))
                          )

    p.chr <- Selected.chr %>% mutate(Method = factor(Method,
                                                     labels=c("Marginal", "Knockoffs"),
                                                     levels=c("Marginal", "Knockoffs"))) %>%
        ggplot(aes(x=BP, y=Importance, color=Method)) +
        geom_point(alpha=1) +
        facet_wrap("Method", ncol=1, scales="free_y") +
        theme_bw() +
        scale_x_continuous(minor_breaks = filter(Selected.knockoffs, CHR==chr)$BP)
        
    # Plot side-by-side
    ggsave(file=sprintf("plots/%s_discoveries_%s_K%d_chr%d.png", phenotype, clumping, K, chr),
           p.chr, width=10, height=7)
}
plot.manhattan(2)
plot.manhattan(6)
plot.manhattan(16)

# Combine results
Results <- Results.marginal %>% filter(Knockoff==FALSE) %>%
    mutate(Importance = -log10(P), Method="Marginal", Sign.W=factor(0, levels=c(-1,0,1))) %>%
    select(-c("Knockoff")) %>%
    full_join(mutate(Results.knockoffs, Method="Knockoffs", Importance=W),
              by=c("CHR", "SNP", "BP", "Group", "Method", "Importance", "Sign.W"))

# Manhattan plot (joint)
threshold.data <- tibble(Method = c("Marginal", "Knockoffs"), Importance = c(-log10(5e-8), W.thres)) %>%
    mutate(Method = factor(Method, levels=c("Marginal", "Knockoffs"), labels=c("Marginal", "Knockoffs")))

p.joint <- Results %>% 
    mutate(Method = factor(Method, levels=c("Marginal", "Knockoffs"),
                           labels=c("Marginal", "Knockoffs")), Sign.W=as.factor(Sign.W)) %>%
    ggplot(aes(x=BP, y=Importance, color=Sign.W)) + geom_point(alpha=0.75) +    
    geom_hline(data=threshold.data, aes(yintercept=Importance), linetype="dashed", color = "black") +
    facet_grid(cols = vars(CHR), rows=vars(Method), scales="free", switch="x") +
    theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0),
          plot.background = element_rect(fill="white"), panel.background = element_rect(fill="gray95"),
          panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(1, "lines"))
ggsave(file=sprintf("plots/%s_manhattan_%s_K%d.png", phenotype, clumping, K), p.joint, width=15, height=7)
