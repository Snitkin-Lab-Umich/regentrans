#test reverse_list_str
z <- list(z1 = list(a = 1, b = 2, c = 3), z2 = list(b = 4, a = 1, c = 0))
o <- list(a = list(z1 = 1, z2 = 1), b = list(z1 = 2, z2 = 4), c = list(z1 = 3, z2 = 0))
oz <- reverse_list_str(z)

test_that("reverse_list_str works", {
  #are o and oz the same?
  #names
  expect_true(all(names(o) == names(oz)))
  #names of z and o are different
  expect_true(length(intersect(names(o), names(z))) == 0)
  #contents of o and oz are the same
  expect_true(all(unlist(oz[[1]]) == unlist(o[[1]])))
  expect_true(all(unlist(oz[[2]]) == unlist(o[[2]])))
  expect_true(all(unlist(oz[[3]]) == unlist(o[[3]])))

})
