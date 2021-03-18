#prune tree function from joyce's github
#goal is to remove second isolate from same patient in same subtree

#first_sts = subtrees(tr)

prune_tree <- function(tree, pt){

  #TODO write checks
  check_prune_tree_input(tree, pt)

  #make sure we reduce the tree to tree tips and pt have in common
  #get the names of the things in common
  isolates <- intersect(tree$tip.label, names(pt))
  #subset locs
  pt_sub <- pt[isolates]
  #subset the tree
  tree_sub <- ape::keep.tip(tree,isolates)
  tree_sub <- phytools::midpoint.root(tree_sub)

  #get subtrees
  first_sts = ape::subtrees(tree_sub)
  #find the order of tree size
  subtr_sizes <- as.numeric(sapply(first_sts, Ntip))
  order(subtr_sizes, decreasing = FALSE)

  #make list of tips to drop
  tips_to_prune <- vector()
  #go into one subtree
  pruned_tree_tips = sapply(1:length(first_sts), FUN = function(x){

    temp_tree = first_sts[[x]]

    #make a list of all of the different isolates
    #for now assume it is facility ID and patient ID
    temp_tree_mat = as.data.frame(sapply(temp_tree$tip.label, FUN = function(x){
      cbind(pt[x])}))
    #if there is more than one isolate from a patient, drop all but the first
    if(!(length(unique(temp_tree_mat[,1])) == nrow(temp_tree_mat))){
      #make a table to show how many per patient
      n_per_pat <- table(temp_tree_mat[,1])
      #subset to ones > 1
      n_per_pat <- n_per_pat[n_per_pat > 1]
      #loop through these patients
      for(i in 1:length(n_per_pat)){
        #find the names of all isolates from that patient
        tips <- rownames(temp_tree_mat)[temp_tree_mat[,1] %in% names(n_per_pat)[i]]
        #find all but first and add them to list to be dropped
        #temp_tree <- drop.tip(temp_tree, tips[2:length(tips)])
        #tips_to_prune <- c(tips_to_prune, tips[2:length(tips)])
        return(tips[2:length(tips)])
      }
    }
  })
  pruned_tips <- unique(unlist(pruned_tree_tips))
  pruned_tree <- drop.tip(tree, pruned_tips)
  return(pruned_tree)
}
