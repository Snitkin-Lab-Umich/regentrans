# get clusters
# tr - tree
# locs - locations #named vector #could be body site, infection vs. colinization, etc. maybe change the name to be groups?
# pureness - how pure the cluster is (<= 0.5)
#bootstrap value should be an argument here



#' Get facility clusters on phylogeny
#'
#' @param tr a tree object returned by the read.tree function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID found in the tree
#' @param pureness how pure each cluster should be (must be > 0.5) (optional, defauly = 1)
#' @param bootstrap Bootstrap support to use to filter unconfident tree edges (optional, default = NULL)
#'
#' @return data.frame of facility clusters on phylogeny
#' @export
#'
#' @examples
get_clusters <- function(tr, locs, pureness = 1, bootstrap = NULL){ # pureness shouldn't be <= 0.5?
  #check inputs
  check_get_clusters_inputs(tr, locs, pureness, bootstrap)

  locs_sub <- locs[tr$tip.label]
  subtrs_sub <- ape::subtrees(tr)
  pure_subtrees <- get_largest_subtree(subtrs = subtrs_sub, isolate_labels = locs_sub, bootstrap = bootstrap, pureness = pureness) # NOTE: this function is in snitkitr right now, but I think we should migrate it to this package (or at least include it here as well); this _might_ be buggy, so definitely good to add unit tests for it
  pure_subtr_info <- bind_cols(f_id=locs_sub,
                               subtr_size=unlist(pure_subtrees$largest_st),
                               index=unlist(pure_subtrees$largest_st_i))
  # change singletons from 0 to 1
  pure_subtr_info <- pure_subtr_info %>% mutate(subtr_size=ifelse(subtr_size==0 & index == 1, 1, subtr_size))
  # remove duplicates (singletons aren't duplicates)
  pure_subtr_info <- pure_subtr_info[!duplicated(pure_subtr_info) | pure_subtr_info$subtr_size == 1,]
  return(pure_subtr_info=pure_subtr_info) #maybe add pureness of cluster?
}

#what does the output mean here? how would I plot it?
