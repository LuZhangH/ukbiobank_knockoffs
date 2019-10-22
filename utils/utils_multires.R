library(gtools)
source("../utils/multilayer_knockoff_filter.R")

prepare_multires <- function(Variants, Stats) {

  ## Cluster structure
  Clusters <- Variants %>% group_by(CHR, Group, Resolution) %>%
    summarise(BP.min=min(BP), BP.max=max(BP)) %>% ungroup() %>%
    select(Resolution, CHR, Group, BP.min, BP.max)

  ## Highest resolution
  res.max <- levels(Variants$Resolution)[length(levels(Variants$Resolution))]

  ## Number of 'variables'
  n <- Clusters %>% filter(Resolution==res.max) %>% nrow()
  ## Number of 'groups'
  M <- length(resolution.list)

  ## Make 'groups'
  groups = matrix(NA, n, M)
  for(m in seq(1,M)) {
    grp.count <- 0
    for(chr in chr.list) {
      Variants.chr <- filter(Variants, CHR==chr, Resolution==resolution.list[m])
      chr.idx <- which(filter(Variants, Resolution==resolution.list[m])$CHR==chr)
      groups[chr.idx,m] <- Variants.chr$Group + grp.count
      grp.count <- grp.count + max(Variants.chr$Group)
    }
  }

  G = apply(groups,2,max) ## G[m] = number of groups for at layer m
  ## for one-way construction, find which groups are "parents" to which other groups
  ## (this construction assumes that for each m>1, groups at layer m
  ## are nested inside those of layer m-1)
  parents <- lapply(2:M, function(m) {
    cat(sprintf("Processing layer %s...\n", m))
    sapply(1:G[m], function(g) {
      idx <- binsearch(function(y) groups[y,m]-g, range=c(1, nrow(groups)))$where[1]
      return(groups[idx,m-1])
    })
  })

  ## Make list of W stats
  W <- lapply(1:M, function(m) {
    cat(sprintf("Processing layer %s...\n", m))
    w.res <- rep(NA, G[m])
    chr.offset <- 0
    for(chr in chr.list) {
      cat(sprintf("Processing chromosome %s...\n", chr))
      chr.idx <- which(filter(Variants, Resolution==resolution.list[m])$CHR==chr)
      Variants.chr <- filter(Variants, CHR==chr, Resolution==resolution.list[m])
      Stats.chr <- filter(Stats, CHR==chr, Resolution==resolution.list[m])
      n.groups <- max(Variants.chr$Group)
      w.chr <- rep(0,n.groups)
      w.chr[Stats.chr$Group] <- Stats.chr$W
      w.res[seq(1,n.groups)+chr.offset] <- w.chr
      chr.offset <- chr.offset + n.groups
    }
    return(w.res)
  })

  ## Return results
  results <- list(groups=groups, parents=parents, W=W)
  return(results)
}
