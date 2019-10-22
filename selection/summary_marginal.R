library(tidyverse)
library(gridExtra)
library(knockoff)

clumping <- "Radj10"
K <- 50
tmp.dir <- "/scratch/PI/candes/ukbiobank_tmp"
phenotype <- "kpolycystic"
regression <- "logistic"

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

Results <- tibble()
for(chr in seq(1,22)) {
    Results.chr <- load.results(chr)
    Results <- rbind(Results, Results.chr)
}

# Compute grouped importance with knockoffs
W.stats <- function(importance, knockoff) {
#    max(importance) * sign(importance[!knockoff]-importance[knockoff])
    importance[!knockoff]-importance[knockoff]
}
Results.knockoffs <- Results %>% group_by(CHR, Group, Knockoff) %>%
    summarise(P=min(P), SNP=SNP[1], BP=round(mean(BP)), Group.size=n()) %>%
    select(CHR, SNP, BP, Knockoff, P, Group, Group.size) %>%
    mutate(Importance = -log10(P)) %>%
    group_by(CHR, Group) %>%
    summarize(W = W.stats(Importance,Knockoff), SNP=SNP[1], BP=BP[1], P=min(P), Group.size=mean(Group.size)) %>%
    mutate(Sign.W=factor(sign(W),levels=c(-1,0,1))) %>%
    arrange(desc(abs(W))) %>%
    select(CHR, SNP, BP, P, Group, Group.size, everything())

# Select groups with knockoffs
W.thres <- knockoff.threshold(Results.knockoffs$W[1:10000], offset=1)
Selected.knockoffs <- Results.knockoffs %>% filter(W > W.thres) %>%
    select(CHR, SNP, BP, Group, P, everything())

# Select groups with marginal testing
Selected.marginal <- Results %>% filter(Knockoff==FALSE) %>%
    group_by(CHR, Group) %>%
    summarise(P=min(P), SNP=SNP[1], BP=round(mean(BP)), Group.size=n()) %>%
    select(CHR, SNP, BP, P, Group, Group.size) %>%
    mutate(Importance = -log10(P)) %>%
    arrange(desc(Importance)) %>%
    filter(P < 5e-8)

# Intersect selections
semi_join(Selected.marginal, Selected.knockoffs, by=c("CHR", "SNP", "BP", "Group"))

# Compare selections on chromosome 16
Selected.knockoffs %>% filter(CHR==16) %>% arrange(P)
Selected.marginal %>% filter(CHR==16) %>% arrange(P)

plot.manhattan <- function(chr) {
    Results.chr <- filter(Results,CHR==chr)
    x.min <- min(Results.chr$BP)
    x.max <- max(Results.chr$BP)

    # Combine discoveries
    Selected.chr <- rbind(Selected.knockoffs %>% filter(CHR==chr) %>%
        mutate(Method = "Knockoffs", Threshold=W.thres, Importance=W),
            Selected.marginal %>% filter(CHR==chr) %>%
        mutate(Method = "Marginal", Threshold=-log10(5e-8))
        )

    p.chr <- Selected.chr %>% mutate(Method = factor(Method, levels=c("Marginal", "Knockoffs"))) %>%
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
plot.manhattan(16)

# Manhattan plot (marginal)
p1 <- Results %>% filter(Knockoff==FALSE) %>%
    mutate(Importance = -log10(P)) %>%
    ggplot(aes(x=BP, y=Importance)) + geom_point(alpha=0.5) +
    geom_hline(yintercept=-log10(5e-8), color="blue") +
    facet_grid(cols = vars(CHR), scales="free_x", switch="x") +
    theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "white", fill = NA, size = 0))

# Manhattan plot (knockoff)
p2 <- Results.knockoffs %>%
    ggplot(aes(x=BP, y=abs(W), color=Sign.W)) + geom_point(alpha=0.5) +
    geom_hline(yintercept=W.thres, color="blue") +
    facet_grid(cols = vars(CHR), scales="free_x", switch="x") +
    theme_bw() +
    theme(legend.position="bottom",
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "white", fill = NA, size = 0))

# Plot side-by-side
p.joint <- arrangeGrob(p1, p2, ncol=1)
ggsave(file=sprintf("plots/%s_manhattan_%s_K%d.png", phenotype, clumping, K), p.joint, width=15, height=7)
