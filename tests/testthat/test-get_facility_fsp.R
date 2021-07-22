#get_facility_fsp output
locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()

test_locs_5 <- locs[locs %in% c("A", "F", "H")]
test_fasta_2 <- aln[names(test_locs_5),]


fsp_mat <- get_facility_fsp(fasta = test_fasta_2, locs = test_locs_5, form = "matrix")
fsp_long <- get_facility_fsp(fasta = test_fasta_2, locs = test_locs_5, form = "long")

test_that("get_facility_fsp works", {
  #for long
  expect_true(class(fsp_long) == "data.frame")
  #ncol = 3
  expect_true(ncol(fsp_long) == 3)
  #colnames
  expect_true(all(colnames(fsp_long) == c("Loc1","Loc2","Fsp")))
  #nrow =
  expect_true(nrow(fsp_long) == (length(unique(test_locs_5))^2)-length(unique(test_locs_5)))
  #col types
  expect_true(all(sapply(fsp_long, class) == c("character","character","numeric")))
  #make sure all vals between 0, 1
  expect_true(all(fsp_long$Fsp <= 1) & all(fsp_long$Fsp > 0))


  #for mat
  expect_true(class(fsp_mat) == "data.frame")
  #nrow = ncol
  expect_true(ncol(fsp_mat) == nrow(fsp_mat))
  #colnames = rownames
  expect_true(all(rownames(fsp_mat) == colnames(fsp_mat)))
  #nrow <= length(locs)
  expect_true(nrow(fsp_mat) <= length(test_locs_5) & nrow(fsp_mat) >= 2)
  #intersect of rownames and locs names is the same group as colnames
  expect_true(length(intersect(rownames(fsp_mat), test_locs_5)) == length(rownames(fsp_mat)))
  #numeric inside
  expect_true(all(fsp_mat <= 1) & all(fsp_mat >= 0))
  #make sure all of the ones on the diagonal are 0
  expect_true(all(diag(as.matrix(fsp_mat)) == 0))

})

