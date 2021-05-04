#test indirect flow
mat <- data.frame(matrix(data = c(0, 20, 12,
                                  20, 0, 26,
                                  30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")
pat_flow <- na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% filter(Var1 != Var2))
colnames(pat_flow) <- c("source_facil", "dest_facil", "n_transfers")
pat_flow$n_transfers <- as.numeric(pat_flow$n_transfers)

test_i_flow_returns <- indirect_flow(pat_flow)
#separate into paths and summary 
test_i_flow_paths <- test_i_flow_returns$paths
test_i_flow_net <- test_i_flow_returns$transfer_network


test_pt_flow_nas <- indirect_flow(pat_flow[-3,])

test_that("indirect_flow works", {
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