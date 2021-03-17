#prune tree function from joyce's github
#goal is to remove second isolate from same patient in same subtree

first_sts = subtrees(tree)

prune tree <- function(tree, pt){

  #TODO write checks

  #make sure we reduce the tree to tree tips and pt have in common
  #get the names of the things in common
  isolates <- intersect(tree$tip.label, names(pt))
  #subset locs
  pt_sub <- pt[isolates]
  #subset the tree
  tree_sub <- ape::keep.tip(tree,isolates)

  #get subtrees
  first_sts = ape::subtrees(tree_sub)

  #go into one subtree
  pruned_tree_tips = sapply(1:length(first_sts), FUN = function(x){print(x)

    temp_tree = first_sts[[x]];

    #make a list of all of the different isolates
    #for now assume it is facility ID and patient ID
    temp_tree_mat = t(sapply(temp_tree$tip.label, FUN = function(x){
      cbind(pt[x],locs[x])}))
    #if there is only one isolate per patient, keep the whole tree
    if(length(unique(temp_tree_mat[,1])) == nrow(temp_tree_mat)){
      keep_tree = temp_tree
    }
    #if there is more than one isolate from a patient, drop all but the first?
    else{
      #make a table to show how many per patient
      n_per_pat <- table(temp_tree_mat[,1])
      #subset to ones > 1
      n_per_pat <- n_per_pat[n_per_pat > 1]
      #loop through these patients
      for(i in 1:length(n_per_pat)){
        #find the indexes of all isolates from that patient
        tips <- which(temp_tree_mat[,1] %in% names(n_per_pat)[i])
        #find all but first
        drop_tips <- tips[2:length(tips)]
        #drop those tips
        temp_tree = drop.tip(temp_tree, rownames(temp_tree_mat)[drop_tips])
      }
      keep_tree <- temp_tree
    }
  }
  #TODO figure out how to get this to the output I want
})

  pruned_tree_tips = unique(unlist(pruned_tree_tips))

  pruned_tree <- drop.tip(tree, pruned_tree_tips)
  return(pruned_tree)
}
