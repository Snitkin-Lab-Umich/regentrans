#' Get facility clusters on the phylogeny
#'
#' @inheritParams get_snv_dists
#' @param tr a tree object returned by the read.tree function from the ape package
#' @param pureness how pure each cluster should be (must be > 0.5) (optional, defauly = 1)
#' @param bootstrap Bootstrap support to use to filter unconfident tree edges (optional, default = NULL)
#'
#' @return list where pure_subtree_info is a data.frame of facility clusters on phylogeny,
#' index indicates which element that cluster is in the list of subtrees,
#' NA indicates it is not part of a subtree;
#' subtrees is an object of the actual subtrees (can be used for plotting);
#' cluster_pureness is the purness of each cluster
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' clusts <- get_clusters(tr, locs, pureness = 1, bootstrap = NULL)
#' }
get_clusters <- function(tr, locs, pureness = 1, bootstrap = NULL){
  #check inputs
  check_get_clusters_inputs(tr, locs, pureness, bootstrap)

  #get the names of the things in common
  isolates <- intersect(tr$tip.label, names(locs))
  #subset locs
  locs_sub <- locs[isolates]
  #subset the tree
  tr <- ape::keep.tip(tr,isolates)

  subtrs_sub <- ape::subtrees(tr)
  pure_subtrees <- get_largest_subtree(subtrs = subtrs_sub, isolate_labels = locs_sub, bootstrap = bootstrap, pureness = pureness)
  pure_subtr_info <- dplyr::bind_cols(loc=locs_sub,
                               subtr_size=unlist(pure_subtrees$largest_st),
                               index=unlist(pure_subtrees$largest_st_i),
                               isolate_id=names(locs_sub))
  # change singletons from 0 to 1
  pure_subtr_info <- pure_subtr_info %>% dplyr::mutate(subtr_size=ifelse(subtr_size==0 & index == 1, 1, subtr_size))
  # change index from 1 to NA
  pure_subtr_info <- pure_subtr_info %>% dplyr::mutate(index=ifelse(index==1, NA, index))
  #add a column to indicate the isolate name if the index = NA
  pure_subtr_info <- pure_subtr_info %>% dplyr::mutate(isolate_id=ifelse(is.na(index), isolate_id, NA))
  # remove duplicates (singletons aren't duplicates)
  pure_subtr_info <- pure_subtr_info[!duplicated(pure_subtr_info$index) | pure_subtr_info$subtr_size == 1,]


  #potential returns if we want to return subtrees too
  returns <- list("pure_subtree_info" = pure_subtr_info, "subtrees" = subtrs_sub, "cluster_pureness" = pureness)

  return(returns)
}


#' Get largest pure subtrees
#'
#' @param subtrs Subtrees created using ape::subtrees to look for clustering on. Should include all isolates of interest.
#' @param isolate_labels Named vector of labels by which pure clusters are defined. Names must be equivalent to tree tip label names.
#' @param control_labels Named vector of labels known to cluster. Names must be equivalent to tree tip label names. This controls for clustering by requiring that the pure clusters must contain multiple of the control labels.
#' @param bootstrap Bootstrap support to use to filter unconfident tree edges (keeps > bootstrap; NULL = keep all; default: 90).
#' @param pureness How pure the subtree has to be to call it a "pure" subtree (default: 1; range 0-1).
#'
#' @return list containing the largest pure subtree that each isolate belongs to, the index of that subtree, and the edges in that subtree.
#' @noRd
#'
get_largest_subtree <- function(subtrs, isolate_labels, control_labels=NULL, bootstrap = 90, pureness = 1){
  #checks
  check_get_largest_subtree_input(subtrs, isolate_labels, control_labels, bootstrap, pureness)

  #largest_st_info = future.apply::future_lapply(names(isolate_labels), function(i){
  largest_st_info = lapply(names(isolate_labels), function(i){
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO.
    #CLUSTERS ARE DEFINED AS:
    # 1) HAVE ONLY A SINGLE EPI LABEL, 2) HAVE BOOTSTRAP SUPPORT GREATER THAN 90, 3) INCLUDE MORE THAN ONE CONTROL LABEL
    #sts = future.apply::future_sapply(subtrs, FUN = function(st){
    sts = sapply(subtrs, FUN = function(st){
      i_in_subtree = i %in% st$tip.label # isolate is in subtree
      st_labs <- isolate_labels[intersect(st$tip.label, names(isolate_labels))]
      one_label = length(unique(st_labs)) == 1 # only one label in subtree
      labs_tab <- table(st_labs)
      labs_max <- max(labs_tab)
      labs_max_lab <- names(labs_tab)[which.max(labs_tab)]
      if(pureness != 1){ # allow some contamination based on pureness threshold
        one_label <- labs_max/length(st_labs) > pureness
      }
      good_bootstrap = rep(TRUE, length(st$node.label[[1]]))
      if(!is.null(bootstrap)){
        good_bootstrap = !is.na(as.numeric(st$node.label[[1]])) && as.numeric(st$node.label[[1]]) > bootstrap # bootstrap support > 90
      }
      multiple_control = ifelse(is.null(control_labels), TRUE, # always true if not controlling for another variable
                                length(unique(control_labels[intersect(st$tip.label, names(control_labels))])) > 1) # more than one control label in subtree
      if(i_in_subtree && one_label && good_bootstrap && multiple_control && isolate_labels[i] == labs_max_lab){
        length(intersect(names(isolate_labels[isolate_labels == labs_max_lab]), st$tip.label))
      }else{
        0
      }
    })
    # GET THE LARGEST SUBTREE
    largest_st = max(sts)
    # GET THE INDEX OF THE LARGEST SUBTREE
    largest_st_i = which.max(sts)
    # GET EDGES BELONGING TO SUBTREES
    largest_st_edges = ape::which.edge(subtrs[[1]], subtrs[[largest_st_i]]$tip.label)
    return(list(largest_st=largest_st,
                largest_st_i=largest_st_i,
                largest_st_edges=largest_st_edges))
  })# end future_apply
  # reverse lists
  largest_st_info = reverse_list_str(largest_st_info)
  return(largest_st_info)
}#end get_largest_subtree


