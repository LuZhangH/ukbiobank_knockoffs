library(Matrix)

############################################################################################
## MULTILAYER KNOCKOFF FILTER

## Arguments
## W:           list of M vectors of knockoff statistics (M is the number of layers)
## groups:      nxM matrix, such that groups[j, m] is the group to which
##              variable j belongs at layer m (n is the number of variables)
## parents:     list of M-1 vectors, such that the m'th vector has g'th entry equal
##              to the index in layer m of the parent of group g at layer m+1
##              (only required for one-way consistency; for two-way consistency may be null)
## q:           length-M vector of FDR target levels
##              (also allowed to be scalar if same across layers)
## offsets:     vector of offsets of length M (also allowed to be scalar if same across layers)
## consistency: the kind of consistency enforced by MKF: either "one-way" or "two-way"
############################################################################################
multilayer_knockoff_filter <- function(W, groups, parents, q, offsets, consistency) {
  ## problem dimensions
  n <- nrow(groups) ## number of variables
  M <- ncol(groups) ## number of layers
  G <- apply(groups,2,max) ## G[m] <- number of groups for at layer m

  cat(sprintf("Running multi-layer filter on %d layers...\n", M))
  
  ## check input for correctness
  stopifnot(length(q) %in% c(1,M))
  stopifnot(length(offsets) %in% c(1, M))
  if(length(q) == 1) {
    q <- rep(q, M)
  }
  if(length(offsets) == 1) {
    offsets <- rep(offsets, M)
  }
  stopifnot(all(apply(groups, 2, function(v)(length(unique(v)))) == G))

  ## reorder groups based on magnitudes of knockoff statistics
  group_orders <- list()               ## ordering of groups at each layer
  groups_reordered <- matrix(0, n, M)  ## reordered group assignments
  parents_reordered <- vector("list", M-1) ## reordered parent assignments
  for(m in 1:M) {
    group_orders[[m]] <- order(abs(W[[m]]), decreasing = TRUE)
    groups_reordered[,m] <- invPerm(group_orders[[m]])[groups[,m]]
    if(m > 1) {
      parents_reordered[[m-1]] <- invPerm(group_orders[[m-1]])[parents[[m-1]][group_orders[[m]]]]
    }
  }

  ## run filter to get threshold indices at each layer
  thresh_idx <- get_thresholds(W, group_orders, groups_reordered, parents_reordered, q, offsets, consistency)

  ## extract knockoff statistic thresholds from threshold indices
  thresh <- numeric(M)
  for(m in 1:M) {
    if(thresh_idx[m] == 0) {
      thresh[m] <- Inf
    }
    else{
      thresh[m] <- abs(W[[m]][group_orders[[m]][thresh_idx[m]]])
    }
  }

  ## define one-bit p-values
  P <- list()
  for(m in 1:M) {
  P[[m]] <- rep(1, G[m])
    kn_stats_ordered <- W[[m]][group_orders[[m]]]
    signs <- sign(kn_stats_ordered)
    P[[m]][signs == 1] <- 0
  }

  ## get selection set
  S_hats_reordered = get_S_hats(P, groups_reordered, parents_reordered, thresh_idx, consistency)
  S_hats = sapply(1:M, function(m)(group_orders[[m]][S_hats_reordered[[m]]]))

  ## return output: selection set, threshold indices, and thresholds
  output <- c()
  output$S_hats <- S_hats
  output$thresh_idx <- thresh_idx
  output$thresh <- thresh
  return(output)
}

## AUXILIARY FUNCTIONS FOR MULTILAYER KNOCKOFF FILTER

## find multilayer knockoff filter thresholds at each layer
get_thresholds <- function(W, group_orders, groups_reordered, parents_reordered, q, offsets, consistency) {
  n <- nrow(groups_reordered)
  M <- ncol(groups_reordered)

  G <- apply(groups_reordered,2,max) ## G[m] <- # groups_reordered, for grouping m

  ## initialize thresholds
  thresh_idx <- G

  ## define one-bit p-values
  P <- list()
  allowable_thresh_idx <- list()
  for(m in 1:M) {
    kn_stats_ordered <- W[[m]][group_orders[[m]]]
    signs <- sign(kn_stats_ordered)
    ## Beginning of old code
    ##allowable_thresh_idx[[m]] <- which(abs(kn_stats_ordered[1:(G[m]-1)]) > abs(kn_stats_ordered[2:G[m]]))
    ## End of old code
    ## Beginning of new code
    if(m==1) {
      sign.flip <- signs[1:(G[m]-1)] > signs[2:G[m]]
      allowable_thresh_idx[[m]] <- which(sign.flip)
    } else {
      #allowable_thresh_idx[[m]] <- 1:G[m] #which(positive)
      allowable_thresh_idx[[m]] = which(abs(kn_stats_ordered[1:(G[m]-1)]) > abs(kn_stats_ordered[2:G[m]]))
    }
    ## End of new code
    allowable_thresh_idx[[m]] <- c(allowable_thresh_idx[[m]], G[m])
    P[[m]] <- rep(1, G[m])
    P[[m]][signs == 1] <- 0
  }

  ## find thresholds
  done <- FALSE
#  while(!done) {
    thresh_idx_old <- thresh_idx
    for(m in 1:M) {
      if(done) break
      if(thresh_idx[m] >= 1) {
        thresh_idx_m <- 0

        cat(sprintf("Computing threshold for layer %d (%d candidate values)...\n",
                    m, length(allowable_thresh_idx[[m]])))
        pb <- txtProgressBar(min = 0, max = length(allowable_thresh_idx[[m]]), style = 3)
        pb.counter <- 0
        for(thresh_idx_m_tmp in rev(allowable_thresh_idx[[m]])) {
          thresh_idx_tmp <- thresh_idx
          thresh_idx_tmp[m] <- thresh_idx_m_tmp
          FDP_hat <- get_FDP_hat(P, groups_reordered, parents_reordered, thresh_idx_tmp, offsets, consistency)
          ##cat(sprintf("Layer %d, %d : %.3f\n", m, thresh_idx_m_tmp, FDP_hat[m]))
          if(FDP_hat[m] <= q[m]) {
            thresh_idx_m <- thresh_idx_m_tmp
            break
          }
          pb.counter <- pb.counter + 1
          setTxtProgressBar(pb, pb.counter)
        }
        close(pb)
        cat(sprintf("Found threshold index %d for layer %d.\n", thresh_idx_m_tmp, m))

        ## Store new threshold
        thresh_idx[m] <- thresh_idx_m

        ## Exit if done
        if(thresh_idx[m] == 0) {
          thresh_idx[m:M] <- 0
          done <- TRUE
          break
        }

        ## Update list of allowed thresholds for the next layer
        if(m<M) {
          ##cat(sprintf("Updating allowable thresholds for layer %d... ", m+1))
          S_hats_tmp <- get_S_hats(P, groups_reordered, parents_reordered, thresh_idx, consistency)
          S_hat_tmp_m <- S_hats_tmp[[m]]
          supported.parents.tmp <- which(parents_reordered[[m]] %in% S_hat_tmp_m)
          allowable_thresh_idx[[m+1]] <- intersect(allowable_thresh_idx[[m+1]], supported.parents.tmp)
          ##cat(sprintf("%d values left.\n", length(allowable_thresh_idx[[m+1]])))
        }
        
      }
      
    }
#    if(all(thresh_idx_old==thresh_idx)) {done <- TRUE}
#  }

  return(thresh_idx)
}

## get estimated number of false discoveries at each layer for given thresholds
get_V_hats <- function(P, thresh_idx, offsets) {
  M <- length(thresh_idx)
  V_hats <- numeric(M)
  V_hats[thresh_idx == 0] <- offsets[thresh_idx == 0]
  for(m in which(thresh_idx > 0)) {
    V_hats[m] <- offsets[m] + sum(P[[m]][1:thresh_idx[m]] == 1)
  }
  return(V_hats)
}

## get selection sets at each layer for given thresholds
get_S_hats <- function(P, groups_reordered, parents_reordered, thresh_idx, consistency) {
  n <- nrow(groups_reordered)
  M <- ncol(groups_reordered)
  if(consistency=="two-way") {
    if(any(thresh_idx == 0)) {
      S_hat <- numeric(0)
    }
    else{
      S_hat <- 1:n ## current selection set
      for(m in 1:M) {
        S_tilde_m <- which(is.element(groups_reordered[,m],intersect(which(P[[m]] == 0), 1:thresh_idx[m])))
        S_hat <- intersect(S_hat, S_tilde_m)
      }
    }
    S_hats <- sapply(1:M, function(m)(unique(groups_reordered[S_hat,m])))
  } else {
    S_hats <- vector("list", M)
    for(m in 1:M) {
      if(thresh_idx[m] == 0) {
        S_hats[[m]] <- numeric(0)
      } else {
        ## Beginning of old code
        ## S_hats[[m]] <- intersect(which(P[[m]] == 0), 1:thresh_idx[m])
        ## End of old code
        ## Beginning of new code
        idx.positive <- which(P[[m]] == 0)
        s.idx <- which(idx.positive <= thresh_idx[m])
        S_hats[[m]] <- idx.positive[s.idx]
        ## End of new code
        if(m > 1) {
          ## Beginning of old code
          ## S_hats[[m]] <- intersect(S_hats[[m]], which(parents_reordered[[m-1]] %in% S_hats[[m-1]]))
          ## End of old code
          ## Beginning of new code
          candidate.parents <- parents_reordered[[m-1]][S_hats[[m]]]
          supported.parents <- which(candidate.parents %in% S_hats[[m-1]])
          S_hats[[m]] <- S_hats[[m]][supported.parents]
          ## End of new code
        }
      }
    }
  }
  return(S_hats)
}

## get FDP-hat for given thresholds
get_FDP_hat <- function(P, groups_reordered, parents_reordered, thresh_idx, offsets, consistency) {
  M <- ncol(groups_reordered)

  S_hats <- get_S_hats(P, groups_reordered, parents_reordered, thresh_idx, consistency)

  V_hats <- numeric(M)
  V_hats[thresh_idx == 0] <- offsets[thresh_idx == 0]
  for(m in which(thresh_idx > 0)) {
    V_hats[m] <- offsets[m] + sum(P[[m]][1:thresh_idx[m]] == 1)
  }

  FDP_hat <- V_hats/sapply(S_hats, function(S_hat)(max(length(S_hat),1)))

  return(FDP_hat)
}
