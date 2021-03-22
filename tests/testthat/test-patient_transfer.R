#tests for patient transfer
#make a source destination pair test matrix
mat <- data.frame(matrix(data = c(0, 20, 30,
                        20, 0, 26,
                        30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")
pat_flow <- na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% filter(Var1 != Var2))
colnames(pat_flow) <- c("source_facil", "dest_facil", "n_transfers")
pat_flow$n_transfers <- as.numeric(pat_flow$n_transfers)

test_locs <- locs[1:4]
test_pt <- as.character(pt[1:4])
names(test_pt) <- names(pt[1:4])
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt)


patient_transfer(pat_flow, test_snv_dists)

