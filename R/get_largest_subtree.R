#' Get largest pure subtrees
#'
#' @param subtrs Subtrees created using ape::subtrees to look for clustering on. Should include all isolates of interest.
#' @param isolate_labels Named vector of labels by which pure clusters are defined. Names must be equivalent to tree tip label names.
#' @param control_labels Named vector of labels known to cluster. Names must be equivalent to tree tip label names. This controls for clustering by requiring that the pure clusters must contain multiple of the control labels.
#' @param bootstrap Bootstrap support to use to filter unconfident tree edges (keeps > bootstrap; NULL = keep all; default: 90).
#' @param pureness How pure the subtree has to be to call it a "pure" subtree (default: 1; range 0-1).
#'
#' @return list containing the largest pure subtree that each isolate belongs to, the index of that subtree, and the edges in that subtree.
#' @export
#'
get_largest_subtree <- function(subtrs, isolate_labels, control_labels=NULL, bootstrap = 90, pureness = 1){
  #checks
  check_get_largest_subtree_input(subtrs, isolate_labels, control_labels, bootstrap, pureness)

  largest_st_info = future.apply::future_lapply(names(isolate_labels), function(i){
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO.
    #CLUSTERS ARE DEFINED AS:
    # 1) HAVE ONLY A SINGLE EPI LABEL, 2) HAVE BOOTSTRAP SUPPORT GREATER THAN 90, 3) INCLUDE MORE THAN ONE CONTROL LABEL
    sts = future.apply::future_sapply(subtrs, FUN = function(st){
      print(i)
      print(st$tip.label)
      i_in_subtree = i %in% st$tip.label # isolate is in subtree
      print(i_in_subtree)
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
