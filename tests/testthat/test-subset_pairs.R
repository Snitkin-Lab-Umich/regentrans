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
