#tests for get_snv_dists
test_locs <- locs[1:3]
test_pt <- as.character(pt[1:3])
names(test_pt) <- names(pt[1:3])
test_dists <- dists[names(test_locs), names(test_locs)]
#one with pt
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt)
#one without pt
test_snv_dist_no_pt <- get_snv_dists(dists = test_dists, locs = test_locs)

mat <- data.frame(matrix(data = c(0, 20, 12,
                                  20, 0, 26,
                                  30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")
test_pt_trans_net <- stats::na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% dplyr::filter(Var1 != Var2))
colnames(test_pt_trans_net) <- c("source_facil", "dest_facil", "n_transfers")
test_pt_trans_net$n_transfers <- as.numeric(test_pt_trans_net$n_transfers)

test_snv_dists_pt_trans <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt, pt_trans_net = test_pt_trans_net)
test_snv_dists_pt_trans_no_pt <- get_snv_dists(dists = test_dists, locs = test_locs, pt = NULL, pt_trans_net = test_pt_trans_net)


test_that("get_snv_dists works", {
  #for with pt
  #check col #s
  expect_true(ncol(test_snv_dists) == 8)
  #check row #s
  expect_true(nrow(test_snv_dists) == (length(test_locs)^2-length(test_locs))/2)
  #check colnames
  expect_true(all(colnames(test_snv_dists) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2", "Patient1", "Patient2", "Pair_Type")))
  #check col types
  expect_true(all(sapply(test_snv_dists, class) == c("factor", "factor", "numeric", "character", "character", "character", "character", "character")))
  # these no longer true when remove duplicate rows (keep only one of s1-s2 s2-s1)
  #check isolates in both lists are in the locs names
  # expect_true(all(setequal(unique(test_snv_dists$Isolate1), unique(test_snv_dists$Isolate2)) & setequal(unique(test_snv_dists$Isolate2), unique(names(test_locs)))))
  #same with patients
  # expect_true(all(setequal(unique(test_snv_dists$Patient1), unique(test_snv_dists$Patient2)) & setequal(unique(test_snv_dists$Patient2), unique(test_pt))))
  #same with locs
  # expect_true(all(setequal(unique(test_snv_dists$Loc1), unique(test_snv_dists$Loc2)) & setequal(unique(test_snv_dists$Loc2), unique(test_locs))))

  #for without pt
  expect_true(ncol(test_snv_dist_no_pt) == 6)
  #check row #s
  expect_true(nrow(test_snv_dist_no_pt) == (length(test_locs)^2-length(test_locs))/2)
  #check colnames
  expect_true(all(colnames(test_snv_dist_no_pt) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2", "Pair_Type")))
  #check col types
  #this won't be true when all of the patients IDs aren't numeric
  expect_true(all(sapply(test_snv_dist_no_pt, class) == c("factor", "factor", "numeric", "character", "character", "character")))
  #check isolates in both lists are in the locs names
  expect_true(all(unique(c(as.character(test_snv_dists$Isolate1), as.character(test_snv_dists$Isolate2))) %in% names(test_locs)))
  # expect_true(all(setequal(unique(test_snv_dists$Isolate1), unique(test_snv_dists$Isolate2)) & setequal(unique(test_snv_dists$Isolate2), unique(names(test_locs)))))
  #same with locs
  expect_true(all(unique(c(as.character(test_snv_dists$Loc1), as.character(test_snv_dists$Loc2))) %in% test_locs))
  # expect_true(all(setequal(unique(test_snv_dists$Loc1), unique(test_snv_dists$Loc2)) & setequal(unique(test_snv_dists$Loc2), unique(test_locs))))

  #one with pt_trans net and pt
  expect_true(ncol(test_snv_dists_pt_trans) == 12)
  #check colnames
  expect_true(all(colnames(test_snv_dists_pt_trans) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2", "Patient1", "Patient2", "Pair_Type", "n_1_to_2_transfers", "n_2_to_1_transfers", "indirect_flow_metric_1_to_2", "indirect_flow_metric_2_to_1")))
  #check coltypes
  expect_true(class(test_snv_dists_pt_trans[,9]) == "numeric")
  expect_true(class(test_snv_dists_pt_trans[,10]) == "numeric")
  #one with pt_trans_net and no pt
  expect_true(ncol(test_snv_dists_pt_trans_no_pt) == 10)

})



#test subset_pairs
test_locs <- locs[1:4]
test_pt <- pt[1:4]
test_dists <- dists[names(test_locs), names(test_locs)]
test_snv_dists_dup <-
  structure(
    list(
      Isolate1 = structure(
        c(2L, 3L, 4L, 1L, 3L, 4L,
          1L, 2L, 4L, 1L, 2L, 3L),
        .Label = c("PCMP_H1", "PCMP_H2", "PCMP_H3",
                   "PCMP_H4"),
        class = "factor"
      ),
      Isolate2 = structure(
        c(1L, 1L,
          1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L),
        .Label = c("PCMP_H1",
                   "PCMP_H2", "PCMP_H3", "PCMP_H4"),
        class = "factor"
      ),
      Pairwise_Dists = c(42,
                         91, 65, 42, 18, 5, 91, 18, 4, 65, 5, 4),
      Loc1 = c("B", "C", "D",
               "A", "C", "D", "A", "B", "D", "A", "B", "C"),
      Loc2 = c("A", "A",
               "A", "B", "B", "B", "C", "C", "C", "D", "D", "D"),
      Patient1 = c(2L,
                   3L, 4L, 1L, 3L, 4L, 1L, 2L, 4L, 1L, 2L, 3L),
      Patient2 = c(1L,
                   1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L),
      Pair_Type = c(
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair",
        "Inter-facility pair"
      )
    ),
    row.names = c(NA,-12L),
    class = "data.frame"
  )

test_subset_pairs <- subset_pairs(test_snv_dists_dup)

test_that("subset_pairs works", {
  #check that there are half the number of isolates
  expect_true(nrow(test_subset_pairs) == nrow(test_snv_dists_dup) * .5)
  #check number of locs vs. rows
  expect_true(nrow(test_subset_pairs) == sum(1:(length(test_locs) - 1)))
  #check that each isolate in locs is represented only length(locs)-1 times
  expect_true(all(table(
    c(test_subset_pairs$Isolate1, test_subset_pairs$Isolate2)
  ) == (length(test_locs) - 1)))
})

