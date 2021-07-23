#tests for get_frac_intra output
locs <- metadata %>% dplyr::select(sample_id, facility) %>% tibble::deframe()

test_locs <- locs[1:3]
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs)
test_snv_df <- structure(list(sample1 = "MN_CRE202", sample2 = "MN_CRE17",
               pairwise_dist = 10, loc1 = "M0211", loc2 = "M0211", pair_type = "intra-facility pair"), row.names = 1958L, class = "data.frame")
test_frac_intra <- get_frac_intra(snv_dists = test_snv_dists)
test_frac_intra_df <- get_frac_intra(snv_dists = test_snv_df)
mat <- data.frame(matrix(data = c(0, 20, 12,
                                  20, 0, 26,
                                  30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")

test_that("get_frac_intra works", {
  #check n cols
  expect_true(ncol(test_frac_intra) == 5)
  #check colnames
  expect_true(all(colnames(test_frac_intra) == c('pairwise_dist','n_intra','n_inter','frac_intra','frac_inter')))
  #check types of all cols
  expect_true(all(sapply(test_frac_intra, class) == c("numeric", "numeric", "numeric", "numeric", "numeric")))
  #check that rows are the divisions of other rows...
  expect_true(all((test_frac_intra$frac_intra == test_frac_intra$n_intra/(test_frac_intra$n_intra + test_frac_intra$n_inter)) | test_frac_intra$frac_intra == 0))
  expect_true(all((test_frac_intra$frac_inter == test_frac_intra$n_inter/(test_frac_intra$n_intra + test_frac_intra$n_inter)) | test_frac_intra$frac_inter == 0))
})
