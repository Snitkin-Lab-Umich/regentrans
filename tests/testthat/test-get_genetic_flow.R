#get_genetic_flow output
locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
pt <- metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()

test_locs_5 <- locs[locs %in% c("A", "F", "H")]
test_pt_5 <- pt[names(pt) %in% names(test_locs_5)]
test_fasta_2 <- aln[names(test_locs_5),]

fsp_mat <- get_genetic_flow(fasta = test_fasta_2, locs = test_locs_5, matrix = TRUE, pt = NULL)
#fsp_pt <- get_genetic_flow(fasta = test_fasta_2, locs = test_locs_5, matrix = TRUE, pt = test_pt_5)
fsp_long <- get_genetic_flow(fasta = test_fasta_2, locs = test_locs_5, matrix = FALSE, pt = NULL)

test_that("get_genetic_flow works", {
  #for long
  expect_true(class(fsp_long) == "data.frame")
  #ncol = 3
  expect_true(ncol(fsp_long) == 3)
  #colnames
  expect_true(all(colnames(fsp_long) == c("loc1","loc2","fsp")))
  #nrow =
  expect_true(nrow(fsp_long) == ((length(unique(test_locs_5))^2)-length(unique(test_locs_5)))/2)
  #col types
  expect_true(all(sapply(fsp_long, class) == c("character","character","numeric")))
  #make sure all vals between 0, 1
  expect_true(all(fsp_long$fsp <= 1) & all(fsp_long$fsp > 0))


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

test_locs_5 <- locs[locs %in% c("A")]
test_fasta_2 <- aln[names(test_locs_5[6:8]),]
test_fasta_2 = test_fasta_2[,apply(test_fasta_2, 2, FUN = function(x){sum(x != x[1] | x == 'n') > 0})]
test_fasta_2 = test_fasta_2[,1:7]
test_fast_char <- as.character(test_fasta_2)
#general test
test_ref <- find_major_alleles(test_fast_char, ref = NULL)

#test if there is a column that is all dashes (it just makes it all dashes)
test_fast_char[,1] = "-"
test_ref_2 <- find_major_alleles(test_fast_char, ref = NULL)

#one with ties
test_fast_char_2 <- rbind(test_fast_char, test_fast_char[3,])
rownames(test_fast_char_2) <- c(rownames(test_fast_char_2)[1:3], "PCMP_H1")
test_ref_3 <- find_major_alleles(test_fast_char_2, ref = NULL)

#one tie with a ref seq
test_ref_4 <- find_major_alleles(test_fast_char_2, ref = c("-", "t", "a", "-", "g", "a", "g"))

test_that("find_major_alleles works", {
  #one normal/easy one with majorities
  expect_true(all(test_ref == c("t", "t", "c", "a", "g", "c", "g")))
  #dashes/missing
  expect_true(all(test_ref_2 == c("-", "t","c", "a", "g", "c", "g")))
  #one tie, should pick the lower-most one
  expect_true(all(test_ref_3 == c("-", "t", "a", "-", "g", "a", "g")))
  #one tie with a ref seq
  expect_true(all(test_ref_4 == c("-", "g", "c", "a", "t", "c", "a")))
})

test_pt_5 <- pt[names(test_locs_5)]
test_fasta_3 <- aln[names(test_pt_5[1:8]),]
test_fasta_3 = test_fasta_3[,apply(test_fasta_3, 2, FUN = function(x){sum(x != x[1] | x == 'n') > 0})]
#test_fast_3_char <- as.character(test_fasta_3[,1:8])
test_meta_seqs <- as.character(make_meta_seqs(test_fasta_3[,1:8], test_locs_5, test_pt_5))

simple_test_meta_seqs <- as.character(make_meta_seqs(test_fasta_3[1:2,1:8], test_locs_5[1:2], test_pt_5[1:2]))

#test where there is a tie
test_fasta_4 <- aln[names(test_pt_5)[2:4],]
test_fasta_4 = test_fasta_4[,apply(test_fasta_4, 2, FUN = function(x){sum(x != x[1] | x == 'n') > 0})]
test_fast_4_char <- as.character(test_fasta_4)
test_meta_seqs_2 <-as.character(make_meta_seqs(test_fasta_4[,50:53], test_locs_5, test_pt_5))

test_that("make_meta_seqs works", {
  #one normal one that doesn't have more than one patient
  expect_true(all(as.character(test_fasta_3[1:2,1:8]) == as.character(simple_test_meta_seqs)))
  #one that has more than one patient
  expect_true(nrow(test_meta_seqs) == length(unique(test_pt_5)))
  #test those sequences
  expect_true((test_meta_seqs[5,] == c("t", "t", "c", "a", "g", "c", "g", "c" )))
  #one that has a tie between the sequences and one is different than the majority alleles
  expect_true(test_meta_seqs_2[1,] == c("c", "a", "t", "a"))

})
