# get clusters
# tr - tree
# locs - locations
# pureness - how pure the cluster is (<= 0.5)
get_clusters <- function(tr, locs, pureness = 1){ # pureness shouldn't be <= 0.5?
  locs_sub <- locs[tr$tip.label]
  subtrs_sub <- subtrees(tr)
  pure_subtrees <- get_largest_subtree(subtrs_sub, locs_sub, bootstrap = NULL, pureness = pureness) # NOTE: this function is in snitkitr right now, but I think we should migrate it to this package (or at least include it here as well); this _might_ be buggy, so definitely good to add unit tests for it
  pure_subtr_info <- bind_cols(f_id=locs_sub,
                               subtr_size=unlist(pure_subtrees$largest_st),
                               index=unlist(pure_subtrees$largest_st_i))
  # change singletons from 0 to 1
  pure_subtr_info <- pure_subtr_info %>% mutate(subtr_size=ifelse(subtr_size==0 & index == 1, 1, subtr_size))
  # remove duplicates (singletons aren't duplicates)
  pure_subtr_info <- pure_subtr_info[!duplicated(pure_subtr_info) | pure_subtr_info$subtr_size == 1,]
  return(pure_subtr_info=pure_subtr_info) #maybe add pureness of cluster?
}
