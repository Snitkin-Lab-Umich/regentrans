# tests for summarize_pairs

locs <- metadata %>% dplyr::select(sample_id, facility) %>% tibble::deframe()
test_locs <- locs[1:3]
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs)
test_isolate_pair_summary <- summarize_pairs(test_snv_dists, summary_fns = c("min"), threshs = c(20, 50))
test_isolate_pair_summary_2 <- summarize_pairs(test_snv_dists, summary_fns = c("min"), threshs = NULL)
test_isolate_pair_summary_3 <- summarize_pairs(test_snv_dists, summary_fns = NULL, threshs = c(20, 50))
test_pt_flow <- pt_trans_df %>% dplyr::filter(source_facil %in% c('C','D','E') & dest_facil %in% c('C','D','E')) %>% get_patient_flow()
test_fsp_long <- make_long_form(fsp)[1:2,]
merged_summaries <- merge_inter_summaries(test_pt_flow, test_isolate_pair_summary, test_fsp_long)
merged_summaries_2 <-  merge_inter_summaries(NULL, test_isolate_pair_summary, test_fsp_long)
merged_summaries_3 <- merge_inter_summaries(test_pt_flow, NULL, test_fsp_long)
merged_summaries_4 <- merge_inter_summaries(test_pt_flow, test_isolate_pair_summary, NULL)

test_that("summarize_pairs works", {
  # everything
  expect_equal(test_isolate_pair_summary,
               structure(list(loc1 = c("A", "A", "B"), loc2 = c("B", "C", "C"
               ), dists_min = c(42, 91, 18), under_20 = c(0L, 0L, 1L), under_50 = c(1L, 0L, 1L)),
               row.names = c(NA, -3L),
               groups = structure(list(loc1 = c("A", "A", "B"), loc2 = c("B", "C", "C"),
                                       .rows = structure(list(1L, 2L, 3L),ptype = integer(0),
                                                         class = c("vctrs_list_of", "vctrs_vctr",
                                                                   "list"))), row.names = c(NA, -3L),
                                  class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
               class = c("grouped_df", "tbl_df", "tbl", "data.frame")))
  # no threshs
  expect_equal(test_isolate_pair_summary_2,
               structure(list(loc1 = c("A", "A", "B"), loc2 = c("B", "C", "C"
               ), dists_min = c(42, 91, 18)), row.names = c(NA, -3L), groups = structure(list(
                 loc1 = c("A", "A", "B"), loc2 = c("B", "C", "C"), .rows = structure(list(
                   1L, 2L, 3L), ptype = integer(0),
                   class = c("vctrs_list_of", "vctrs_vctr", "list"))), row.names = c(NA, -3L),
                 class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
               class = c("grouped_df", "tbl_df", "tbl", "data.frame")))
  # no summary_fns
  expect_equal(test_isolate_pair_summary_3,
               structure(list(loc1 = c("A", "A", "B"), loc2 = c("B", "C", "C"),
                              under_20 = c(0L, 0L, 1L), under_50 = c(1L, 0L, 1L)),
               class = "data.frame", row.names = c(NA, -3L)))
})

test_that("merge_inter_summaries works", {
  expect_equal(merged_summaries,
               structure(list(loc1 = c("A", "A", "B", "A"),
                              loc2 = c("B", "C", "C", "D"),
                              dists_min = c(42, 91, 18, NA),
                              under_20 = c(0L, 0L, 1L, NA),
                              under_50 = c(1L, 0L, 1L, NA),
                              fsp = c(NA, 0.767963355449144,
                                      NA, 0.746959511800763),
                              n_transfers_f12 = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              pt_trans_metric_f12 = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              n_transfers_f21 = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              pt_trans_metric_f21 = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              sum_transfers = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              sum_pt_trans_metric = c(NA_real_, NA_real_, NA_real_, NA_real_)),
                         row.names = c(NA, -4L),
                         groups = structure(list(loc1 = c("A", "A", "A", "B"),
                                                 loc2 = c("B", "C", "D", "C"),
                                                 .rows = structure(list(1L, 2L, 4L, 3L),
                                                                   ptype = integer(0),
                                                                   class = c("vctrs_list_of",
                                                                             "vctrs_vctr", "list"))),
                                            row.names = c(NA, -4L), class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
                         class = c("grouped_df","tbl_df", "tbl", "data.frame")))
  expect_equal(merged_summaries_2,
               structure(list(loc1 = c("A", "A", "B", "A"),
                              loc2 = c("B", "C", "C", "D"),
                              dists_min = c(42, 91, 18, NA),
                              under_20 = c(0L, 0L, 1L, NA),
                              under_50 = c(1L, 0L, 1L, NA),
                              fsp = c(NA, 0.767963355449144,
                                      NA, 0.746959511800763)),
                         row.names = c(NA, -4L),
                         groups = structure(list(loc1 = c("A", "A", "A", "B"),
                                                 loc2 = c("B", "C", "D", "C"),
                                                 .rows = structure(list(1L, 2L, 4L, 3L),
                                                                   ptype = integer(0),
                                                                   class = c("vctrs_list_of",
                                                                             "vctrs_vctr", "list"))),
                                            row.names = c(NA, -4L),
                                            class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
                         class = c("grouped_df", "tbl_df", "tbl", "data.frame")))
  expect_equal(merged_summaries_3,
               structure(list(loc1 = c("A", "A"), loc2 = c("C", "D"),
                              fsp = c(0.767963355449144,
                                      0.746959511800763),
                              n_transfers_f12 = c(NA_real_, NA_real_),
                              pt_trans_metric_f12 = c(NA_real_, NA_real_),
                              n_transfers_f21 = c(NA_real_, NA_real_),
                              pt_trans_metric_f21 = c(NA_real_, NA_real_),
                              sum_transfers = c(NA_real_, NA_real_),
                              sum_pt_trans_metric = c(NA_real_, NA_real_)),
                         class = "data.frame", row.names = c(NA, -2L)))
  expect_equal(merged_summaries_4,
               structure(list(loc1 = c("A", "A", "B"), loc2 = c("B", "C", "C"), dists_min = c(42, 91, 18), under_20 = c(0L, 0L, 1L), under_50 = c(1L, 0L, 1L), n_transfers_f12 = c(NA_real_, NA_real_, NA_real_), pt_trans_metric_f12 = c(NA_real_, NA_real_, NA_real_), n_transfers_f21 = c(NA_real_, NA_real_, NA_real_), pt_trans_metric_f21 = c(NA_real_, NA_real_, NA_real_), sum_transfers = c(NA_real_, NA_real_, NA_real_), sum_pt_trans_metric = c(NA_real_, NA_real_, NA_real_)), row.names = c(NA, -3L), groups = structure(list(loc1 = c("A", "A", "B"), loc2 = c("B", "C", "C"), .rows = structure(list(1L, 2L, 3L), ptype = integer(0), class = c("vctrs_list_of", "vctrs_vctr", "list"))), row.names = c(NA, -3L), class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE), class = c("grouped_df",  "tbl_df", "tbl", "data.frame")))
})
