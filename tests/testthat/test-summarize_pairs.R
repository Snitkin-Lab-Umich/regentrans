# tests for summarize_pairs

locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
pt <- metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()
test_locs <- locs[locs %in% c("C", "D", "E")]
test_pt <- pt[names(pt) %in% names(test_locs)]
test_dists <- dists[names(test_locs), names(test_locs)]
test_pair_types <- get_pair_types(dists = test_dists, locs = test_locs, pt = test_pt)
test_isolate_pair_summary <- summarize_pairs(test_pair_types, summary_fns = c("min"), threshs = c(20, 50))
test_isolate_pair_summary_2 <- summarize_pairs(test_pair_types, summary_fns = c("min"), threshs = NULL)
test_isolate_pair_summary_3 <- summarize_pairs(test_pair_types, summary_fns = NULL, threshs = c(20, 50))
test_pt_flow <- regentrans::pt_trans_df %>% dplyr::filter(source_facil %in% c('C','D','E') & dest_facil %in% c('C','D','E')) %>% get_patient_flow()
test_fsp_long <- make_long_form(fsp) %>% dplyr::filter(loc1 %in% c('C','D','E'), loc2 %in% c('C','D','E'))
merged_summaries <- merge_inter_summaries(test_pt_flow, test_isolate_pair_summary, test_fsp_long)
merged_summaries[is.na(merged_summaries)] <- 0
merged_summaries <- merged_summaries %>% dplyr::mutate(dplyr::across(is.numeric, round, digits = 2))
merged_summaries_2 <-  merge_inter_summaries(NULL, test_isolate_pair_summary, test_fsp_long)
merged_summaries_2[is.na(merged_summaries_2)] <- 0
merged_summaries_2 <- merged_summaries_2 %>% dplyr::mutate(dplyr::across(is.numeric, round, digits = 2))
merged_summaries_3 <- merge_inter_summaries(test_pt_flow, NULL, test_fsp_long)
merged_summaries_3 <- merged_summaries_3 %>% dplyr::mutate(dplyr::across(is.numeric, round, digits = 2))
merged_summaries_4 <- merge_inter_summaries(test_pt_flow, test_isolate_pair_summary, NULL)
merged_summaries_4[is.na(merged_summaries_4)] <- 0
merged_summaries_4 <- merged_summaries_4 %>% dplyr::mutate(dplyr::across(is.numeric, round, digits = 2))


test_that("summarize_pairs works", {
  # everything
  expect_equal(test_isolate_pair_summary,
               structure(list(loc1 = c("C", "C", "C", "D", "D", "E"),
                              loc2 = c("C", "D", "E", "D", "E", "E"),
                              dist_min = c(0, 4, 2, 0, 6, 0),
                              leq_20 = c(16, 44, 40, 110, 29, 95),
                              leq_50 = c(101, 319, 175, 327, 251, 232)),
               row.names = c(NA, -6L),
               groups = structure(list(loc1 = c("C", "C", "C", "D", "D", "E"), loc2 = c("C", "D", "E", "D", "E", "E"),
                                       .rows = structure(list(1L, 2L, 3L, 4L, 5L, 6L),ptype = integer(0),
                                                         class = c("vctrs_list_of", "vctrs_vctr",
                                                                   "list"))), row.names = c(NA, -6L),
                                  class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
               class = c("grouped_df", "tbl_df", "tbl", "data.frame")))
  # no threshs
  expect_equal(test_isolate_pair_summary_2,
               structure(list(loc1 = c("C", "C", "C", "D", "D", "E"),
                              loc2 = c("C", "D", "E", "D", "E", "E"),
                              dist_min = c(0, 4, 2, 0, 6, 0)), row.names = c(NA, -6L),
                              groups = structure(list(
                 loc1 = c("C", "C", "C", "D", "D", "E"), loc2 = c("C", "D", "E", "D", "E", "E"), .rows = structure(list(
                   1L, 2L, 3L, 4L, 5L, 6L), ptype = integer(0),
                   class = c("vctrs_list_of", "vctrs_vctr", "list"))), row.names = c(NA, -6L),
                 class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
               class = c("grouped_df", "tbl_df", "tbl", "data.frame")))
  # no summary_fns
  expect_equal(test_isolate_pair_summary_3,
               structure(list(loc1 = c("C", "C", "C", "D", "D", "E"),
                              loc2 = c("C", "D", "E", "D", "E", "E"),
                              leq_20 = c(16, 44, 40, 110, 29, 95),
                              leq_50 = c(101, 319, 175, 327, 251, 232)),
               class = "data.frame", row.names = c(NA, -6L)))
})

test_that("merge_inter_summaries works", {
  expect_true(all(merged_summaries ==
               structure(list(loc1 = c("C", "C", "C","C", "C", "D", "D", "D", "E"),
                              loc2 = c("C", "D", "D", "E", "E", "D", "E", "E", "E"),
                              dist_min = c(0, 4, 4, 2, 2, 0, 6, 6, 0),
                              leq_20 = c(16, 44, 44, 40, 40, 110, 29, 29, 95),
                              leq_50 = c(101, 319, 319, 175, 175, 327, 251, 251, 232),
                              fsp = c(0, 0.34, 0.34, 0.34, 0.34, 0, 0.24, 0.24, 0),
                              n_transfers_f12 = c(0, 474, 474, 1718, 1718, 0, 504, 504, 0),
                              pt_trans_metric_f12 = c(0, 0.22, 0.22, 0.78, 0.78, 0, 0.66, 0.66, 0),
                              n_transfers_f21 = c(0, 2660, 2660, 1632, 1632, 0, 330,  330,   0),
                              pt_trans_metric_f21 = c(0, 0.84, 0.84, 0.83, 0.83, 0, 0.18, 0.18, 0),
                              sum_transfers = c(0, 3134, 3134, 3350, 3350, 0,  834, 834, 0),
                              sum_pt_trans_metric = c(0, 1.06, 1.06, 1.62, 1.62, 0, 0.84, 0.84, 0)),
                         row.names = c(NA, -9L),
                         groups = structure(list(loc1 = c("C", "C", "C","C", "C", "D", "D", "D", "E"),
                                                 loc2 = c("C", "D", "D", "E", "E", "D", "E", "E", "E"),
                                                 .rows = structure(list(1L, 2L, 4L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
                                                                   ptype = integer(0),
                                                                   class = c("vctrs_list_of",
                                                                             "vctrs_vctr", "list"))),
                                            row.names = c(NA, -9L), class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
                         class = c("grouped_df","tbl_df", "tbl", "data.frame"))))
  expect_true(all(merged_summaries_2 ==
               structure(list(loc1 = c("C", "C", "C", "D", "D", "E"),
                              loc2 = c("C", "D", "E", "D", "E", "E"),
                              dist_min = c(0, 4, 2, 0, 6, 0),
                              leq_20 = c(16, 44, 40, 110, 29, 95),
                              leq_50 = c(101, 319, 175, 327, 251, 232),
                              fsp = c(0, 0.34, 0.34, 0, 0.24, 0)),
                         row.names = c(NA, -6L),
                         groups = structure(list(loc1 = c("C", "C", "C", "D", "D", "E"),
                                                 loc2 = c("C", "D", "E", "D", "E", "E"),
                                                 .rows = structure(list(1L, 2L, 4L, 3L, 4L, 5L, 6L),
                                                                   ptype = integer(0),
                                                                   class = c("vctrs_list_of",
                                                                             "vctrs_vctr", "list"))),
                                            row.names = c(NA, -6L),
                                            class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
                         class = c("grouped_df", "tbl_df", "tbl", "data.frame"))))
  expect_true(all(merged_summaries_3 ==
               structure(list(loc1 = c( "C", "C", "C", "C", "D", "D"),
                              loc2 = c("D", "D", "E", "E", "E", "E"),
                              fsp = c(0.34, 0.34, 0.34, 0.34, 0.24, 0.24),
                              n_transfers_f12 = c(474,  474, 1718, 1718,  504,  504),
                              pt_trans_metric_f12 = c(0.22, 0.22, 0.78, 0.78, 0.66, 0.66),
                              n_transfers_f21 = c(2660, 2660, 1632, 1632,  330,  330),
                              pt_trans_metric_f21 = c(0.84, 0.84, 0.83, 0.83, 0.18, 0.18),
                              sum_transfers = c(3134, 3134, 3350, 3350,  834,  834),
                              sum_pt_trans_metric = c(1.06, 1.06, 1.62, 1.62, 0.84, 0.84)),
                         class = "data.frame", row.names = c(NA, -6L))))
  expect_true(all(merged_summaries_4 ==
               structure(list(loc1 = c("C", "C", "C","C", "C", "D", "D", "D", "E"),
                              loc2 = c("C", "D", "D", "E", "E", "D", "E", "E", "E"),
                              dist_min = c(0, 4, 4, 2, 2, 0, 6, 6, 0),
                              leq_20 = c(16, 44, 44, 40, 40, 110, 29, 29, 95),
                              leq_50 = c(101, 319, 319, 175, 175, 327, 251, 251, 232),
                              n_transfers_f12 = c(0, 474, 474, 1718, 1718, 0, 504, 504, 0),
                              pt_trans_metric_f12 = c(0, 0.22, 0.22, 0.78, 0.78, 0, 0.66, 0.66, 0),
                              n_transfers_f21 = c(0, 2660, 2660, 1632, 1632, 0, 330,  330,   0),
                              pt_trans_metric_f21 = c(0, 0.84, 0.84, 0.83, 0.83, 0, 0.18, 0.18, 0),
                              sum_transfers = c(0, 3134, 3134, 3350, 3350, 0,  834, 834, 0),
                              sum_pt_trans_metric = c(0, 1.06, 1.06, 1.62, 1.62, 0, 0.84, 0.84, 0)),
                         row.names = c(NA, -9L),
                         groups = structure(list(loc1 = c("C", "C", "C","C", "C", "D", "D", "D", "E"),
                                                 loc2 = c("C", "D", "D", "E", "E", "D", "E", "E", "E"),
                                                 .rows = structure(list(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
                                                                   ptype = integer(0),
                                                                   class = c("vctrs_list_of", "vctrs_vctr", "list"))),
                                            row.names = c(NA, -9L), class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
                         class = c("grouped_df",  "tbl_df", "tbl", "data.frame"))))
})
