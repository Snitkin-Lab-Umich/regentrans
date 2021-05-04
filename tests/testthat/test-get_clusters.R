#test get_clusters output
test_locs <- regentrans::locs[1:10]
test_tr <- ape::keep.tip(regentrans::tr,names(test_locs))
test_clusters <- get_clusters(tr = test_tr,locs = test_locs)
test_pure_subtree_info <- test_clusters$pure_subtree_info
test_subtrees <- test_clusters$subtrees
test_cluster_pureness <- test_clusters$cluster_pureness
test_first_subtr <- test_subtrees[1][[1]]
st_tiplabs <- sapply(1:nrow(test_clusters$pure_subtree_info), function(x){
  i <- test_clusters$pure_subtree_info$index[x]
  name <- test_clusters$pure_subtree_info$isolate_name[x]
  if(!is.na(i)){
    name <- test_clusters$subtrees[[i]]$tip.label
  }
  name
})

test_that("get_snv_dists works", {
  #check that test_clusters is a list of three
  expect_true(length(test_clusters) == 3)
  #check the types
  expect_true(any(class(test_pure_subtree_info) == "data.frame"))
  expect_true(class(test_cluster_pureness) == "numeric")
  expect_true(class(test_subtrees) == "list")
  #make sure it is a list of trees (types of elements in the list are phylo)
  expect_true(all(unname(sapply(test_subtrees, class)) == "phylo"))
  #check test_pure_subtree info ncol
  expect_true(ncol(test_pure_subtree_info) == 4)
  #colnames
  expect_true(all(colnames(test_pure_subtree_info) == c("f_id", "subtr_size", "index", "isolate_name")))
  #nrow, at most one for each isolate, at least one (if they are all in the same subtree)
  expect_true((nrow(test_pure_subtree_info) <= length(test_locs)) & (nrow(test_pure_subtree_info) >= 1))
  #check column types
  expect_true(all(sapply(test_pure_subtree_info, class) == c("character", "numeric", "integer", "character")))
  #make sure if there is an index there is no isolate name
  expect_true(all(filter(test_pure_subtree_info, !is.na(test_pure_subtree_info$index))$isolate_name == " "))
  #make sure if there is.na(index) size == 1 and isolate name != " "
  expect_true(all(filter(test_pure_subtree_info, is.na(test_pure_subtree_info$index))$isolate_name != " ") &
                all(filter(test_pure_subtree_info, is.na(test_pure_subtree_info$index))$subtr_size == 1))
  #make sure all of the fids exist in the test_locs
  #this might break down when we add tree pruning??
  expect_true(length(intersect(unique(test_locs), unique(test_pure_subtree_info$f_id))) ==
                length(unique(test_pure_subtree_info$f_id)) &
                length(intersect(unique(test_locs), unique(test_pure_subtree_info$f_id))) ==
                length(unique(test_locs)))
  #make sure first subtree is the whole tree
  expect_true(all(test_first_subtr$tip.label == test_tr$tip.label))
  #make sure the sum of the subtree lengths is the length of the isolate list
  expect_true(sum(test_pure_subtree_info$subtr_size) == length(test_locs))
  #make sure the tip labels match
  expect_true(length(unlist(st_tiplabs)) == ape::Ntip(test_tr))
})
