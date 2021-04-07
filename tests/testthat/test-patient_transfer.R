#tests for patient transfer
#make a source destination pair test matrix
mat <- data.frame(matrix(data = c(0, 20, 12,
                        20, 0, 26,
                        30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")
pat_flow <- na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% filter(Var1 != Var2))
colnames(pat_flow) <- c("source_facil", "dest_facil", "n_transfers")
pat_flow$n_transfers <- as.numeric(pat_flow$n_transfers)

test_locs <- locs[1:3]
test_pt <- as.character(pt[1:3])
names(test_pt) <- names(pt[1:3])
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt)

test_pt_trans <- patient_transfer(pat_flow, test_snv_dists, thresh = 50)

test_that("patient_transfer works", {
  #check ncols
  expect_true(ncol(test_pt_trans) == 5)
  #check colnames
  expect_true(all(colnames(test_pt_trans) == c("Loc1", "Loc2", "n_closely_related_pairs", "n_1_to_2_transfers", "n_2_to_1_transfers")))
  #check coltypes
  expect_true(all(sapply(test_pt_trans, class) == c("character", "character", "integer", "numeric", "numeric")))
  #check that there are <= nrow pat_dlow
  expect_true(nrow(test_pt_trans) <= nrow(pat_flow))
})
