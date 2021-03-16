#test subset_pairs 
test_locs <- locs[1:4]
test_pt <- pt[1:4]
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt)

test_subset_pairs <- subset_pairs(test_snv_dists)

test_that("subset_pairs works", {
  #check that there are half the number of isolates
  expect_true(nrow(test_subset_pairs) == nrow(test_snv_dists)*.5)
  #check number of locs vs. rows 
  expect_true(nrow(test_subset_pairs) == sum(1:(length(test_locs)-1)))
  #check that each isolate in locs is represented only length(locs)-1 times 
  expect_true(all(table(c(test_subset_pairs$Isolate1, test_subset_pairs$Isolate2)) == (length(test_locs)-1)))
})