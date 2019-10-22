#!/usr/bin/env Rscript
.libPaths("/home/groups/candes/Software/miniconda2/envs/ukb/lib/R/library")

library(tidyverse)

## Parameters
chr <- 22
resolution <- "Radj100"
seed <- 0

scratch <- "/scratch/PI/candes/ukbiobank_tmp"

## Load MAF table
maf.file <- sprintf("%s/knockoff_diagnostics/%s_K50_s%s/ukb_gen_chr%s.frq", scratch, resolution, seed, chr)
Maf <- read_table(maf.file, col_types=cols())

## Load knockoff key table
key.file <- sprintf("%s/knockoffs/%s_K50_s%s/ukb_gen_chr%s.key", scratch, resolution, seed, chr)
Key <- read_delim(key.file, delim=" ", col_types=cols())

## Cross-reference key with MAF
Maf.original <- Maf %>%
  inner_join(transmute(Key, CHR=Chr, SNP=Variant, BP=Position, Knockoff), by = c("CHR", "SNP")) %>%
  filter(Knockoff==F) %>%
  mutate(SNP=str_replace(SNP, ".A", ""), SNP=str_replace(SNP, ".B", "")) %>%
  select(CHR,SNP,BP,MAF)
Maf.knockoffs <- Maf %>%
  inner_join(transmute(Key, CHR=Chr, SNP=Variant, BP=Position, Knockoff), by = c("CHR", "SNP")) %>%
  filter(Knockoff==T) %>%
  mutate(SNP=str_replace(SNP, ".A", ""), SNP=str_replace(SNP, ".B", "")) %>%
  select(CHR,SNP,BP,MAF)
Maf.diagnostics <- Maf.original %>% inner_join(Maf.knockoffs, by=c("CHR", "SNP", "BP"))
colnames(Maf.diagnostics) <- c("CHR", "SNP", "BP", "MAF.x", "MAF.xk")

## Load LD table
ld.file <- sprintf("%s/knockoff_diagnostics/%s_K50_s%s/ukb_gen_chr%s.ld", scratch, resolution, seed, chr)
LD <- read_table(ld.file, col_types=cols())

## Load group information
grp.file <- sprintf("%s/clumping/%s/grp_chr%d.txt", scratch, resolution, chr)
Clumping <- read_delim(grp.file, delim=" ", col_types=cols())

## Add group and knockoff info to LD table
Key.tmp <- Key %>% transmute(CHR_A=Chr, SNP_A=Variant, Group_A=Group, Knockoff_A=Knockoff)
LD <- LD %>% left_join(Key.tmp, by = c("CHR_A", "SNP_A"))
Key.tmp <- Key %>% transmute(CHR_B=Chr, SNP_B=Variant, Group_B=Group, Knockoff_B=Knockoff)
LD <- LD %>% left_join(Key.tmp, by = c("CHR_B", "SNP_B"))
LD <- LD %>% mutate(SNP_A=str_replace(SNP_A, ".A", ""), SNP_A=str_replace(SNP_A, ".B", ""),
                    SNP_B=str_replace(SNP_B, ".A", ""), SNP_B=str_replace(SNP_B, ".B", ""))
LD <- LD %>% select(SNP_A, SNP_B, Group_A, Group_B, Knockoff_A, Knockoff_B, BP_A, BP_B, R2)

## Compute LD diagnostics
Stats <- LD %>% filter(Group_A!=Group_B) %>%
  group_by(SNP_A, SNP_B, BP_A, BP_B) %>%
  summarise(R2.xk=mean(R2[which(Knockoff_A*Knockoff_B==T)]),
            R2.x=mean(R2[which((!Knockoff_A)*(!Knockoff_B)==T)]),
            R2.xxk=mean(R2[which((Knockoff_A)*(!Knockoff_B)==T|(!Knockoff_A)*(Knockoff_B)==T)])
            ) %>%
  ungroup() %>%
  mutate(BP_diff=abs(BP_A-BP_B))
Stats$R2.x[is.nan(Stats$R2.x)] <- 0
Stats$R2.xk[is.nan(Stats$R2.xk)] <- 0
Stats$R2.xxk[is.nan(Stats$R2.xxk)] <- 0

## Compute LD diagnostics (self)
Stats.self <- LD %>% filter(SNP_A==SNP_B) %>%
  select(SNP_A, SNP_B, BP_A, BP_B, R2)

## Add MAF info to LD stats
Stats.maf <- Stats %>%
  left_join(transmute(Maf.original, CHR=CHR, SNP_A=SNP, MAF=MAF)) %>%
  left_join(transmute(Maf.original, CHR=CHR, SNP_B=SNP, MAF=MAF))

## Save results
out.maf.file <- sprintf("%s/plots/gof_maf_chr%d_%s.txt", scratch, chr, resolution)
Maf.diagnostics %>% write_delim(out.maf.file, delim=" ")
cat(sprintf("MAF diagnostics written to: %s\n", out.maf.file))

out.ld.file <- sprintf("%s/plots/gof_ld_chr%d_%s.txt", scratch, chr, resolution)
Stats.maf %>% write_delim(out.ld.file, delim=" ")
cat(sprintf("LD diagnostics written to: %s\n", out.ld.file))

out.self.file <- sprintf("%s/plots/gof_self_chr%d_%s.txt", scratch, chr, resolution)
Stats.self %>% write_delim(out.self.file, delim=" ")
cat(sprintf("LD (self) diagnostics written to: %s\n", out.self.file))

if(FALSE) {

  ## Plot MAF
  Maf.diagnostics %>%
    ggplot(aes(x=MAF.x, y=MAF.xk)) +
    geom_point(alpha=0.1) +
    xlim(0,0.5) + ylim(0,0.5) +
    theme_bw()

  ## Plot LD table
  Stats %>%
    ggplot(aes(x=R2.x, y=R2.xk)) +
    geom_point(alpha=0.1) +
    xlim(0,1) + ylim(0,1) +
    theme_bw()

  Stats %>%
    ggplot(aes(x=R2.x, y=R2.xxk)) +
    geom_point(alpha=0.1) +
    xlim(0,1) + ylim(0,1) +
    theme_bw()

  Stats.maf %>%
    filter(MAF>=0.05) %>%
    ggplot(aes(x=R2.x, y=R2.xk)) +
    geom_point(alpha=0.1) +
    xlim(0,1) + ylim(0,1) +
    theme_bw()

  Stats.maf %>%
    filter(MAF>=0.05) %>%
    ggplot(aes(x=R2.x, y=R2.xxk)) +
    geom_point(alpha=0.1) +
    xlim(0,1) + ylim(0,1) +
    theme_bw()
}
