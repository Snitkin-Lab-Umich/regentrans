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
#' @return list where pure_subtree_info is a data.frame of facility clusters on phylogeny, index indicates which element that cluster is in the list of subtrees, NA indicates it is not part of a subtree; subtrees is an object of the actual subtrees (can be used for plotting);cluster_pureness is the purness of each cluster
#' @export
#'
#' @examples
get_clusters <- function(tr, locs, pureness = 1, bootstrap = NULL){ # pureness shouldn't be <= 0.5?
  #check inputs
  check_get_clusters_inputs(tr, locs, pureness, bootstrap)

  #get the names of the things in common
  isolates <- intersect(tr$tip.label, names(locs))

  locs_sub <- locs[isolates]
  subtrs_sub <- ape::subtrees(tr)
  pure_subtrees <- get_largest_subtree(subtrs = subtrs_sub, isolate_labels = locs_sub, bootstrap = bootstrap, pureness = pureness) # NOTE: this function is in snitkitr right now, but I think we should migrate it to this package (or at least include it here as well); this _might_ be buggy, so definitely good to add unit tests for it
  pure_subtr_info <- bind_cols(f_id=locs_sub,
                               subtr_size=unlist(pure_subtrees$largest_st),
                               index=unlist(pure_subtrees$largest_st_i),
                               isolate_name=names(locs_sub))
  # change singletons from 0 to 1
  pure_subtr_info <- pure_subtr_info %>% mutate(subtr_size=ifelse(subtr_size==0 & index == 1, 1, subtr_size))
  # change index from 1 to NA
  pure_subtr_info <- pure_subtr_info %>% mutate(index=ifelse(index==1, NA, index))
  #add a column to indicate the isolate name if the index = NA
  pure_subtr_info <- pure_subtr_info %>% mutate(isolate_name=ifelse(is.na(index), isolate_name, " "))
  # remove duplicates (singletons aren't duplicates)
  pure_subtr_info <- pure_subtr_info[!duplicated(pure_subtr_info$index) | pure_subtr_info$subtr_size == 1,]


  #potential returns if we want to return subtrees too
  returns <- list("pure_subtree_info" = pure_subtr_info, "subtrees" = pure_subtrees, "cluster_pureness" = pureness)

  return(returns)
}


