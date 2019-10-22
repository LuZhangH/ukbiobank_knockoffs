loci.intersect <- function(bp.min.x, bp.max.x, bp.min.y, bp.max.y, gap=1e5) {
    intersect <- as.logical((bp.min.y<=bp.max.x+gap) * (bp.max.y>=bp.min.x-gap))
    return(intersect)
}

find_supported <- function(Lower.in, Upper.in) {
    Lower <- Lower.in %>% select(Resolution, CHR, BP.min, BP.max)
    Upper <- Upper.in %>% select(Resolution, CHR, BP.min, BP.max)
    Supported <- full_join(Lower,Upper,by="CHR") %>%
        filter(BP.min.y>=BP.min.x, BP.max.y<=BP.max.x) %>%
        transmute(Resolution=Resolution.y, CHR=CHR, BP.min=BP.min.y, BP.max=BP.max.y)
    Supported <- inner_join(Supported, Upper.in, by = c("Resolution", "CHR", "BP.min", "BP.max"))
    return(Supported)
}

find_support <- function(Lower.in, Upper.in) {
    Lower <- Lower.in %>% select(Resolution, CHR, BP.min, BP.max)
    Upper <- Upper.in %>% select(Resolution, CHR, BP.min, BP.max)
    Supported <- full_join(Lower,Upper,by="CHR") %>%
        filter(BP.min.y>=BP.min.x, BP.max.y<=BP.max.x) %>%
        transmute(Resolution=Resolution.x, CHR=CHR, BP.min=BP.min.x, BP.max=BP.max.x)
    Supported <- inner_join(Supported, Lower.in, by = c("Resolution", "CHR", "BP.min", "BP.max"))
    return(Supported)
}

remove_top <- function(Lower.in, Upper.in) {
    Lower <- Lower.in %>% select(Resolution, CHR, BP.min, BP.max)
    Upper <- Upper.in %>% select(Resolution, CHR, BP.min, BP.max)
    Supported <- full_join(Lower,Upper,by="CHR") %>%
        filter(BP.min.y>=BP.min.x, BP.max.y<=BP.max.x) %>%
        transmute(Resolution=Resolution.x, CHR=CHR, BP.min=BP.min.x, BP.max=BP.max.x)
    Outer <- anti_join(Lower.in, Supported, by = c("Resolution", "CHR", "BP.min", "BP.max"))
    return(Outer)
}

consistent_filter <- function(Selections, resolution.min=2) {
    if(nrow(Selections)==0){
        return(Selections)
    }
    resolution.list <- sort(unique(Selections$Resolution))
    Selections.supported <- Selections %>% filter(Resolution==resolution.min)
    if(nrow(Selections.supported)==0) {
        return(Selections[0,])
    }
    k.list <- 1:pmax(1,length(resolution.list)-1) 
    for(k in k.list) {
        resolution <- resolution.list[k]
        resolution.up <- resolution.list[k+1]
        Lower <- Selections.supported %>% filter(Resolution==resolution)
        Upper <- Selections %>% filter(Resolution==resolution.up)
        Upper.supported <- find_supported(Lower, Upper)
        Selections.supported <- rbind(Selections.supported, Upper.supported)
    }
    Selections.supported <- Selections.supported %>% distinct(Resolution, CHR, Group, .keep_all = TRUE)
    return(Selections.supported)
}

consistent_to_outer <- function(Selections) {
    if(nrow(Selections)==0){
        return(Selections)
    }
    resolution.list <- sort(unique(Selections$Resolution))
    Selections.outer <- Selections %>% filter(Resolution==resolution.list[length(resolution.list)])
    if(length(resolution.list)==1) {
        return(Selections.outer)
    }
    for(k in seq(length(resolution.list)-1,1)) {
        resolution <- resolution.list[k]
        resolution.up <- resolution.list[k+1]
        Lower <- Selections %>% filter(Resolution==resolution)
        Upper <- Selections %>% filter(Resolution==resolution.up)
        Lower.outer <- remove_top(Lower, Upper)
        Selections.outer <- rbind(Selections.outer, Lower.outer)
    }
    Selections.outer <- Selections.outer %>% distinct(Resolution, CHR, Group, .keep_all = TRUE)
    return(Selections.outer)
}

outer_to_consistent <- function(Outer, Selections) {    
    resolution.list <- sort(unique(Selections$Resolution))
    # List findings at the highest resolution
    Consistent <- Outer %>% filter(Resolution==max(Outer$Resolution))
    # List findings at the lower resolutions
    for(k in order(resolution.list,decreasing=T)) {
        resolution <- resolution.list[k]
        if(resolution < max(Outer$Resolution) ) {
            resolution.up <- resolution.list[k+1]
            # Find selections below consistent
            Upper <- Consistent %>% filter(Resolution==resolution.up)
            Lower <- Selections %>% filter(Resolution==resolution)
            Lower.support <- find_support(Lower, Upper)            
            Lower.outer <- Outer %>% filter(Resolution==resolution)
            Consistent <- rbind(Consistent, Lower.support)
            Consistent <- rbind(Consistent, Lower.outer)
        }
    }
    Consistent <- Consistent %>% distinct(Resolution, CHR, Group, .keep_all = TRUE)
    return(Consistent)
}
