#tests for get_frac_intra output
test_locs <- locs[1:3]
test_pt <- as.character(pt[1:3])
names(test_pt) <- names(pt[1:3])
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt)
test_snv_df <- structure(list(Isolate1 = "MN_CRE202", Isolate2 = "MN_CRE17",
               Pairwise_Dists = 10, Loc1 = "M0211", Loc2 = "M0211", Pair_Type = "Intra-facility pair"), row.names = 1958L, class = "data.frame")
test_frac_intra <- get_frac_intra(snv_dists = test_snv_dists)
test_frac_intra_df <- get_frac_intra(snv_dists = test_snv_df)
mat <- data.frame(matrix(data = c(0, 20, 12,
                                  20, 0, 26,
                                  30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")
test_pt_trans_net <- stats::na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% dplyr::filter(Var1 != Var2))
colnames(test_pt_trans_net) <- c("source_facil", "dest_facil", "n_transfers")
test_pt_trans_net$n_transfers <- as.numeric(test_pt_trans_net$n_transfers)
test_snv_dists_pt_trans_net <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt, pt_trans_net = test_pt_trans_net)

test_frac_intra_pt_trans_net <- get_frac_intra(snv_dists = test_snv_dists_pt_trans_net)

test_that("get_frac_intra works", {
  #check n cols
  expect_true(ncol(test_frac_intra) == 5)
  #check colnames
  expect_true(all(colnames(test_frac_intra) == c('Pairwise_Dists','n_Intra','n_Inter','Frac_Intra','Frac_Inter')))
  #check types of all cols
  expect_true(all(sapply(test_frac_intra, class) == c("numeric", "numeric", "numeric", "numeric", "numeric")))
  #check that rows are the divisions of other rows...
  expect_true(all((test_frac_intra$Frac_Intra == test_frac_intra$n_Intra/(test_frac_intra$n_Intra + test_frac_intra$n_Inter)) | test_frac_intra$Frac_Intra == 0))
  expect_true(all((test_frac_intra$Frac_Inter == test_frac_intra$n_Inter/(test_frac_intra$n_Intra + test_frac_intra$n_Inter)) | test_frac_intra$Frac_Inter == 0))

  #check it works with the input with pt trans net
  #check n cols
  expect_true(ncol(test_frac_intra_pt_trans_net) == 5)
  #check colnames
  expect_true(all(colnames(test_frac_intra_pt_trans_net) == c('Pairwise_Dists','n_Intra','n_Inter','Frac_Intra','Frac_Inter')))
  #check types of all cols
  expect_true(all(sapply(test_frac_intra_pt_trans_net, class) == c("numeric", "numeric", "numeric", "numeric", "numeric")))
  #check that rows are the divisions of other rows...
  expect_true(all((test_frac_intra_pt_trans_net$Frac_Intra == test_frac_intra_pt_trans_net$n_Intra/(test_frac_intra_pt_trans_net$n_Intra + test_frac_intra_pt_trans_net$n_Inter)) | test_frac_intra_pt_trans_net$Frac_Intra == 0))
  expect_true(all((test_frac_intra_pt_trans_net$Frac_Inter == test_frac_intra_pt_trans_net$n_Inter/(test_frac_intra_pt_trans_net$n_Intra + test_frac_intra_pt_trans_net$n_Inter)) | test_frac_intra_pt_trans_net$Frac_Inter == 0))
  # check that get_snv_dists part works
  expect_equal(expect_message(get_frac_intra(dists = test_dists, locs = test_locs, pt = test_pt),
                              'Running get_snv_dists...'),
    test_frac_intra)
})
