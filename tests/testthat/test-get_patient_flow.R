#tests for patient transfer
#make a source destination pair test matrix
locs <- metadata %>% dplyr::select(sample_id, facility) %>% tibble::deframe()

mat <- data.frame(matrix(data = c(0, 20, 12,
                        20, 0, 26,
                        30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")
pat_flow <- stats::na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% dplyr::filter(Var1 != Var2))
colnames(pat_flow) <- c("source_facil", "dest_facil", "n_transfers")
pat_flow$n_transfers <- as.numeric(pat_flow$n_transfers)
pt_flow_sub <- pat_flow[1:3,]

test_locs <- locs[1:3]
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs)
#one without paths returned
test_pt_trans <- get_patient_flow(edge_df = pat_flow)
#one with paths returned
test_pt_trans_paths <- get_patient_flow(pat_flow, paths = TRUE)

test_that("get_patient_flow works", {
  #for without paths return
  #check ncols
  expect_true(ncol(test_pt_trans) == 8)
  #check colnames
  expect_true(all(colnames(test_pt_trans) == c("loc1", "loc2", "n_transfers_f12", "pt_trans_metric_f12", "n_transfers_f21",
                                               "pt_trans_metric_f21", "sum_transfers", "sum_pt_trans_metric")))
  #check coltypes
  expect_true(all(sapply(test_pt_trans, class) == c(loc1 = "character", loc2 = "character", n_transfers_f12 = "numeric",
                                                    pt_trans_metric_f12 = "numeric", n_transfers_f21 = "numeric",
                                                    pt_trans_metric_f21 = "numeric", sum_transfers = "numeric", sum_pt_trans_metric = "numeric"
  )))
  #check that there are <= nrow pat_dlow
  expect_true(nrow(test_pt_trans) <= nrow(pat_flow))

  #check paths return
  #check that it is a list
  expect_true(class(test_pt_trans_paths) == "list")
  #length 2
  expect_true(length(test_pt_trans_paths) == 2)
  #list names
  expect_true(all(names(test_pt_trans_paths) == c("pt_trans_summary", "paths")))
  #list types
  expect_true(all(sapply(test_pt_trans_paths, class) == c("data.frame", "list")))
  # duplicate source/destination facility rows returns error
  expect_error(get_patient_flow(dplyr::bind_rows(pat_flow,pat_flow),test_snv_dists),
               "Multiple rows in the patient transfer network contain the same source and destination facility. Please include only unique source and destination pairs.")
  # missing patient transfers for some facilities
  expect_warning(expect_equal(get_patient_flow(pt_flow_sub),
               structure(list(loc1 = c("A", "A", "B"), loc2 = c("B", "C", "C"),
               n_transfers_f12 = c(NA_real_, NA_real_, NA_real_),
               pt_trans_metric_f12 = c(0, 0, 0), n_transfers_f21 = c(20, 12, NA),
               pt_trans_metric_f21 = c(1, 1, 0), sum_transfers = c(NA_real_, NA_real_, NA_real_),
               sum_pt_trans_metric = c(1, 1, 0)), row.names = c(NA, -3L), class = "data.frame")))
  })


#test indirect flow
mat <- data.frame(matrix(data = c(0, 20, 12,
                                  20, 0, 26,
                                  30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")
pat_flow <- stats::na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% dplyr::filter(Var1 != Var2))
colnames(pat_flow) <- c("source_facil", "dest_facil", "n_transfers")
pat_flow$n_transfers <- as.numeric(pat_flow$n_transfers)

test_i_flow_returns <- get_indirect_flow(pat_flow)
#separate into paths and summary
test_i_flow_paths <- test_i_flow_returns$paths
test_i_flow_net <- test_i_flow_returns$transfer_network


test_pt_flow_nas <- get_indirect_flow(pat_flow[-3,])

test_that("get_indirect_flow works", {
  #test a normal one (whole return)
  expect_true(class(test_i_flow_returns) == "list")
  expect_true(length(test_i_flow_returns) == 2)
  expect_true(all(names(test_i_flow_returns) == c("transfer_network", "paths")))
  #list types
  expect_true(any(class(test_i_flow_returns$transfer_network) == "data.frame"))
  expect_true(class(test_i_flow_returns$paths) == "list")
  #test normal one (paths)
  #length
  expect_true(length(test_i_flow_returns$paths) == 4)
  #names
  expect_true(all(names(test_i_flow_returns$paths) == c("vpath", "epath", "predecessors", "inbound_edges")))

  #test normal one (network)
  #class
  expect_true(any(class(test_i_flow_net) == "data.frame"))
  #ncol
  expect_true(ncol(test_i_flow_net) == 3)
  #colnames
  expect_true(all(colnames(test_i_flow_net) == c("source_facil", "dest_facil", "pt_trans_metric")))
  #coltypes
  expect_true(all(sapply(test_i_flow_net, class) == c("character", "character", "numeric")))
  #nrow
  expect_true(nrow(test_i_flow_net) == nrow(pat_flow))

  #test one with some NAs...
  expect_true(sum(is.na(test_pt_flow_nas$transfer_network$pt_trans_metric)) == 0)


})

