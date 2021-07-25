#test get_clusters output
locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()

test_locs <- locs[1:10]
test_tr <- ape::keep.tip(tr,names(test_locs))
test_clusters <- get_clusters(tr = test_tr,locs = test_locs)
test_pure_subtree_info <- test_clusters$pure_subtree_info
test_subtrees <- test_clusters$subtrees
test_cluster_pureness <- test_clusters$cluster_pureness
test_first_subtr <- test_subtrees[1][[1]]
st_tiplabs <- sapply(1:nrow(test_clusters$pure_subtree_info), function(x){
  i <- test_clusters$pure_subtree_info$index[x]
  name <- test_clusters$pure_subtree_info$isolate_id[x]
  if(!is.na(i)){
    name <- test_clusters$subtrees[[i]]$tip.label
  }
  name
})

test_that("get_pair_types works", {
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
  expect_true(all(colnames(test_pure_subtree_info) == c("loc", "subtr_size", "index", "isolate_id")))
  #nrow, at most one for each isolate, at least one (if they are all in the same subtree)
  expect_true((nrow(test_pure_subtree_info) <= length(test_locs)) & (nrow(test_pure_subtree_info) >= 1))
  #check column types
  expect_true(all(sapply(test_pure_subtree_info, class) == c("character", "numeric", "integer", "character")))
  #make sure if there is an index there is no isolate name
  expect_true(all(is.na(dplyr::filter(test_pure_subtree_info, !is.na(test_pure_subtree_info$index))$isolate_id)))
  #make sure if there is.na(index) size == 1 and isolate name != " "
  expect_true(all(dplyr::filter(test_pure_subtree_info, is.na(test_pure_subtree_info$index))$isolate_id != " ") &
                all(dplyr::filter(test_pure_subtree_info, is.na(test_pure_subtree_info$index))$subtr_size == 1))
  #make sure all of the fids exist in the test_locs
  #this might break down when we add tree pruning??
  expect_true(length(intersect(unique(test_locs), unique(test_pure_subtree_info$loc))) ==
                length(unique(test_pure_subtree_info$loc)) &
                length(intersect(unique(test_locs), unique(test_pure_subtree_info$loc))) ==
                length(unique(test_locs)))
  #make sure first subtree is the whole tree
  expect_true(all(test_first_subtr$tip.label == test_tr$tip.label))
  #make sure the sum of the subtree lengths is the length of the isolate list
  expect_true(sum(test_pure_subtree_info$subtr_size) == length(test_locs))
  #make sure the tip labels match
  expect_true(length(unlist(st_tiplabs)) == ape::Ntip(test_tr))
})



#tests for get_largest_subtree
test_locs <- locs[1:12]
test_tr <- ape::keep.tip(tr,names(test_locs))
test_subtr <- ape::subtrees(test_tr)
test_pure_subtrees <- get_largest_subtree(subtrs = test_subtr, isolate_labels = test_locs, bootstrap = NULL, pureness = 1)
test_pure_subtrees_2 <- get_largest_subtree(subtrs = test_subtr, isolate_labels = test_locs, bootstrap = NULL, pureness = 0.6)
test_pure_subtrees_3 <- get_largest_subtree(subtrs = test_subtr, isolate_labels = test_locs, bootstrap = 100, pureness = 1)
p1_i <- unlist(test_pure_subtrees$largest_st_i)
p2_i <- unlist(test_pure_subtrees_2$largest_st_i)
p3_i <- unlist(test_pure_subtrees_3$largest_st_i)


test_that("get_largest_subtree works", {
  #test that it is a list of length three
  expect_true(length(test_pure_subtrees) == 3)
  #test the names of items in the list
  expect_true(all(names(test_pure_subtrees) == c("largest_st","largest_st_i","largest_st_edges")))
  #test the types of the items in the list
  expect_true(all(sapply(test_pure_subtrees, class) == c("list", "list", "list")))
  expect_true(all(sapply(test_pure_subtrees[[1]], class) == "numeric"))
  expect_true(all(sapply(test_pure_subtrees[[2]], class) == "integer"))
  expect_true(all(sapply(test_pure_subtrees[[3]], class) == "integer"))
  #test that the number of isolates is the same between locs and in the list
  expect_true(all(sapply(test_pure_subtrees, length) == length(test_locs)))
  #test for isolates we largest subtree != 0 that their index != 1
  expect_true(all((test_pure_subtrees[[1]] != 0) == (test_pure_subtrees[[2]] != 1)))
  expect_equal(p1_i, c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 3L, 3L, 1L, 1L, 1L))
  expect_equal(p2_i, c(1L, 1L, 7L, 1L, 1L, 1L, 1L, 3L, 3L, 1L, 1L, 7L))
  expect_equal(p3_i, c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L))
})

