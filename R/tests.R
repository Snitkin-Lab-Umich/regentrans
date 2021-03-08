#tests for regentrans - to be made when the public facing dataset is available
####################################test data######################################
#test_locs <- locs[1:4]
#test_pt <- pt[1:4]
#test_dists <- dists[names(test_locs), names(test_locs)]


####################################test checks######################################
# get_snv_dists
test_that("check_get_snv_dists_input_pt works", {
  #check ones that should pass
  #with pt
  expect_true(is.null(check_get_snv_dists_input_pt(test_dists, test_locs, test_pt)))
  #without pt
  expect_true()

  expect_error(
    check_dataset("not_a_df"),
    "The dataset must be a `data.frame` or `tibble`"
  )
  expect_error(
    check_dataset(data.frame(outcome = c(), var1 = c())),
    "No rows detected in dataset."
  )
  expect_error(
    check_dataset(data.frame(outcome = 1:3)),
    "1 or fewer columns detected in dataset. There should be an outcome column and at least one feature column."
  )
})

##################################test functions#####################################
