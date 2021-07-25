#tests for get_pair_types
locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
test_locs <- locs[1:3]
test_dists <- dists[names(test_locs), names(test_locs)]
test_pair_types <- get_pair_types(dists = test_dists, locs = test_locs)


mat <- data.frame(matrix(data = c(0, 20, 12,
                                  20, 0, 26,
                                  30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")

test_that("get_pair_types works", {
  # check number of columns
  expect_true(ncol(test_pair_types) == 6)
  #check row #s
  expect_true(nrow(test_pair_types) == (length(test_locs)^2-length(test_locs))/2)
  #check colnames
  expect_true(all(colnames(test_pair_types) == c("isolate1", "isolate2", "pairwise_dist", "loc1", "loc2", "pair_type")))
  #check col types
  #this won't be true when all of the patients IDs aren't numeric
  expect_true(all(sapply(test_pair_types, class) == c("factor", "factor", "numeric", "character", "character", "character")))
  #check isolates in both lists are in the locs names
  expect_true(all(unique(c(as.character(test_pair_types$isolate1), as.character(test_pair_types$isolate2))) %in% names(test_locs)))
  # expect_true(all(setequal(unique(test_pair_types$isolate1), unique(test_pair_types$isolate2)) & setequal(unique(test_pair_types$isolate2), unique(names(test_locs)))))
  #same with locs
  expect_true(all(unique(c(as.character(test_pair_types$loc1), as.character(test_pair_types$loc2))) %in% test_locs))
  # expect_true(all(setequal(unique(test_pair_types$loc1), unique(test_pair_types$loc2)) & setequal(unique(test_pair_types$loc2), unique(test_locs))))

})



#test subset_pairs
test_locs <- locs[1:4]
test_dists <- dists[names(test_locs), names(test_locs)]
test_pair_types_dup <-
  structure(
    list(
      isolate1 = structure(
        c(2L, 3L, 4L, 1L, 3L, 4L,
          1L, 2L, 4L, 1L, 2L, 3L),
        .Label = c("PCMP_H1", "PCMP_H2", "PCMP_H3",
                   "PCMP_H4"),
        class = "factor"
      ),
      isolate2 = structure(
        c(1L, 1L,
          1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L),
        .Label = c("PCMP_H1",
                   "PCMP_H2", "PCMP_H3", "PCMP_H4"),
        class = "factor"
      ),
      pairwise_dist = c(42,
                         91, 65, 42, 18, 5, 91, 18, 4, 65, 5, 4),
      loc1 = c("B", "C", "D",
               "A", "C", "D", "A", "B", "D", "A", "B", "C"),
      loc2 = c("A", "A",
               "A", "B", "B", "B", "C", "C", "C", "D", "D", "D"),
      Patient1 = c(2L,
                   3L, 4L, 1L, 3L, 4L, 1L, 2L, 4L, 1L, 2L, 3L),
      Patient2 = c(1L,
                   1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L),
      pair_type = c(
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair",
        "inter-facility pair"
      )
    ),
    row.names = c(NA,-12L),
    class = "data.frame"
  )

test_subset_pairs <- subset_pairs(test_pair_types_dup)

test_that("subset_pairs works", {
  #check that there are half the number of isolates
  expect_true(nrow(test_subset_pairs) == nrow(test_pair_types_dup) * .5)
  #check number of locs vs. rows
  expect_true(nrow(test_subset_pairs) == sum(1:(length(test_locs) - 1)))
  #check that each isolate in locs is represented only length(locs)-1 times
  expect_true(all(table(
    c(test_subset_pairs$isolate1, test_subset_pairs$isolate2)
  ) == (length(test_locs) - 1)))
})

