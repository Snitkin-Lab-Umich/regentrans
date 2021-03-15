#tests for get_frac_intra output 
test_locs <- locs[1:4]
test_pt <- as.character(pt[1:4])
names(test_pt) <- names(pt[1:4])
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt)
test_threshs <- seq(1,50,1)
test_frac_intra <- get_frac_intra(snv_dists = test_snv_dists, threshs = test_threshs)


test_that("get_frac_intra works", {
  #check n cols
  expect_true(ncol(test_frac_intra) == 5)
  #check that n rows matches threshs 
  expect_true(nrow(test_frac_intra) <= length(test_threshs))
  #check colnames 
  expect_true(all(colnames(test_frac_intra) == c('Thresh','n_Intra','n_Inter','Frac_Intra','Frac_Inter')))
  #check types of all cols 
  expect_true(all(sapply(test_frac_intra, class) == c("numeric", "numeric", "numeric", "numeric", "numeric")))
  #check that rows are the divisions of other rows... 
  expect_true(all(test_frac_intra$Frac_Intra == test_frac_intra$n_Intra/(test_frac_intra$n_Intra + test_frac_intra$n_Inter)))
  expect_true(all(test_frac_intra$Frac_inter == test_frac_intra$n_Inter/(test_frac_intra$n_Intra + test_frac_intra$n_Inter)))
  #check that the length of the overlap between threshs and the column is less than or equal nrow 
  expext_true(length(intersect(test_frac_intra$Thresh, test_threshs)) == nrow(test_frac_intra))
})