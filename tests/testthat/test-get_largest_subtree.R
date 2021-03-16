#tests for get_largest_subtree
test_locs <- locs[1:10]
test_tr <- ape::keep.tip(tr,names(test_locs))
test_subtr <- ape::subtrees(test_tr)
test_pure_subtrees <- get_largest_subtree(subtrs = test_subtr, isolate_labels = test_locs, bootstrap = NULL, pureness = 1)

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
})

