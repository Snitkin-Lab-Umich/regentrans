###################################test data######################################
locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
pt <- metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()

test_locs <- locs[1:4]
test_locs_2 <- locs[5:8]
test_locs_3 <- locs[1:3]
test_locs_4 <- locs[1]
test_locs_5 <- locs[locs %in% c("A", "F", "H")]
test_locs_6 <- locs[1:4]
test_locs_7 <- test_locs_5[1:4]
test_locs_8 <- test_locs_5[1:5]
test_locs_9 <- test_locs_5[1:6]
test_locs_fsp <- test_locs
test_locs_fsp[1:4] <- "A"
test_pt <- pt[1:4]
test_pt_2 <- pt[5:8]
test_pt_3 <- pt[1:3]
test_pt_5 <- pt[names(pt) %in% names(test_locs_5)]
test_pt_6 <- pt[names(pt) %in% names(test_locs_6)]
test_pt_7 <- pt[names(pt) %in% names(test_locs_7)]
test_pt_8 <- pt[names(pt) %in% names(test_locs_8)]
test_pt_9 <- pt[names(pt) %in% names(test_locs_9)]
test_dists <- dists[names(test_locs), names(test_locs)]
test_dists_2 <- dists[names(test_locs_3), names(test_locs_3)]
test_dists_3 <- test_dists_2
test_dists_3[2] <- -1
test_dists_4 <- test_dists_2
test_dists_4[,1] <- as.character(test_dists_4[,1])
test_pair_types <- get_pair_types(dists = test_dists, locs = test_locs, pt = test_pt)
test_pair_types_2 <- test_pair_types[,2:ncol(test_pair_types)]
test_pair_types_3 <- test_pair_types
colnames(test_pair_types_3) <- c("A", "B", "C", "D", "E", "F")
test_pair_types_4 <- test_pair_types
test_pair_types_4$pairwise_dist <- as.character(test_pair_types_4$pairwise_dist)
test_pair_types_5 <- test_pair_types %>% tibble::tibble() %>% dplyr::filter(isolate1 == 'a') %>% data.frame()
test_tr <- ape::keep.tip(tr,names(test_locs))
test_tr_2 <- ape::keep.tip(tr,names(test_locs_3))
test_fasta <- aln[names(test_locs),]
test_fasta_2 <- aln[names(test_locs_5),]
test_fasta_3 <- aln[names(test_locs_8),]
ls <- list(list(a = 2, b = 3), list(c = "a", d = "b"))
ls_2 <- list(list(a = 2, b = 3), list(c = "a", d = "b"), "x")
test_subtr <- ape::subtrees(test_tr)
test_subtr_2 <- ape::subtrees(test_tr_2)
test_pt_flow <- regentrans::pt_trans_df %>% dplyr::filter(source_facil %in% c('C','D','E') & dest_facil %in% c('C','D','E')) %>% get_patient_flow()
test_isolate_pair_summary <- summarize_pairs(test_pair_types)
test_fsp_long <- make_long_form(fsp)

mat <- data.frame(matrix(data = c(0, 20, 30, 40,
                                  20, 0, 26, 50,
                                  30, 26, 0, 57,
                                  40, 50, 57, 0), nrow = 4, ncol = 4))
rownames(mat) <- c("A", "B", "C", "D")
colnames(mat) <- c("A", "B", "C", "D")
test_pt_trans_net <- stats::na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% dplyr::filter(Var1 != Var2))
colnames(test_pt_trans_net) <- c("source_facil", "dest_facil", "n_transfers")
test_pt_trans_net$n_transfers <- as.numeric(test_pt_trans_net$n_transfers)
test_pt_trans_net_2 <- test_pt_trans_net[,2:ncol(test_pt_trans_net)]
test_pt_trans_net_3 <- test_pt_trans_net
colnames(test_pt_trans_net_3) <- c("A", "B", "C")
test_pt_trans_net_4 <- test_pt_trans_net
test_pt_trans_net_4$n_transfers <- as.character(test_pt_trans_net_4$n_transfers)
test_pt_trans_net_5 <- test_pt_trans_net %>% dplyr::filter(source_facil == "A", dest_facil == "A")

test_pair_types_pt_trans <- get_pair_types(dists = test_dists, locs = test_locs, pt = test_pt)

mat2 <- data.frame(matrix(data = c(0, 20, 12,
                                  20, 0, 26,
                                  30, 26, 0), nrow = 3, ncol = 3))
rownames(mat2) <- c("A", "B", "C")
colnames(mat2) <- c("A", "B", "C")
test_pt_trans_net2 <- stats::na.omit(data.frame(as.table(as.matrix(mat2))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% dplyr::filter(Var1 != Var2))
colnames(test_pt_trans_net2) <- c("source_facil", "dest_facil", "n_transfers")
test_pt_trans_net2$n_transfers <- as.numeric(test_pt_trans_net2$n_transfers)

mat3 <- mat
rownames(mat3) <- c(1, 2, 3, 4)

####################################test get_pair_types######################################
test_that("check_get_pair_types_input works", {
  #check ones that should pass
  expect_true(is.null(check_get_pair_types_input(test_dists, test_locs, test_pt)))
  #one input - dists
  expect_error(
    check_get_pair_types_input(dists = test_dists),
    'argument "locs" is missing, with no default',
    fixed = TRUE
  )
  #one input - locs
  expect_error(
    check_get_pair_types_input(locs = test_locs),
    'argument "dists" is missing, with no default',
    fixed = TRUE
  )
  #no pt input
  expect_error(
    check_get_pair_types_input(dists = test_dists, locs = test_locs),
    'argument "pt" is missing, with no default',
    fixed = TRUE
  )
  #test_dists not a dists object
  expect_error(
    check_get_pair_types_input(dists = "test_dists", locs = test_locs, pt = test_pt),
    "The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided: character",
    fixed = TRUE
  )
  #locs not a named list
  expect_error(
    check_get_pair_types_input(dists = test_dists, locs = "test_locs", pt = test_pt),
    "The locs object must be a named list of locations named by sample IDs",
    fixed = TRUE
  )
  #dist names not matching locs names
  expect_error(
    check_get_pair_types_input(dists = test_dists, locs = test_locs_2, pt = test_pt),
    "You have provided an isolate patient ID vector that does not contain all isolates found in the patient location vector. Please subset",
    fixed = TRUE
  )
  #one with less than two in common
  expect_error(
    check_get_pair_types_input(dists = test_dists, locs = test_locs_4, pt = test_pt),
    "You have only supplied locations for 1 isolates. Please supply a named vector of locations for at least 2 isolates",
    fixed = TRUE
  )
  #one where pt is null
  expect_true(is.null(check_get_pair_types_input(test_dists, test_locs, NULL)))
  #one where pt is the wrong type
  expect_error(
    check_get_pair_types_input(dists = test_dists, locs = test_locs, pt = "test_pt"),
    "The pt object must be a named list of patient IDs named by sample IDs or NULL",
    fixed = TRUE
  )
  #one where pt has fewer than locs
  expect_error(
    check_get_pair_types_input(dists = test_dists, locs = test_locs, pt = test_pt[1:3]),
    "You have supplied a vector of more isolates (n =  4 ) with locations than exist in your patient ID vector (n =  3 ). Please subset",
    fixed = TRUE
  )
  #one where pt has different than locs
  expect_error(
    check_get_pair_types_input(dists = test_dists, locs = test_locs, pt = test_pt_2),
    "You have provided an isolate patient ID vector that does not contain all isolates found in the patient location vector. Please subset",
    fixed = TRUE
  )

  #warnings
  #dist names being a subset of locs names
  warnings <- capture_warnings(check_get_pair_types_input(dists = test_dists_2, locs = test_locs, pt = test_pt))
  expect_true(warnings[1] == "You have supplied a list of more isolates (n =  4 ) with locations than exist in your SNV distance matrix (n =  3 ). Will subset")
  expect_true(warnings[2] == "You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset")
  #locs names being a subset of dists names
  warnings <- capture_warnings(check_get_pair_types_input(dists = test_dists, locs = test_locs_3, pt = test_pt_3))
  expect_true(warnings[1] == "You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset")
  warnings <- capture_warnings(check_get_pair_types_input(dists = test_dists_2, locs = test_locs, pt = test_pt))
  expect_true(warnings[1] == "You have supplied a list of more isolates (n =  4 ) with locations than exist in your SNV distance matrix (n =  3 ). Will subset")
  expect_true(warnings[2] == "You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset")

})

# check_dists

test_that("check_dists works", {
  expect_error(check_dists(dists = test_dists_2[,1:2]),
               "The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but the dimensions of your matrix are not equal. The matrix you provided is 3 x 2")
  expect_error(check_dists(
    dists = data.frame(test_dists_2[1,1])),
    "Your SNV matrix only has  1  samples. Please use a SNV distance matrix that includes 2 or more samples")
  expect_error(check_dists(
    dists = (test_dists_2[1,1])),
    "The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided: numeric")
  expect_error(check_dists(
    dists = test_dists_3),
    "The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided an object with values < 0")
  expect_error(check_dists(test_dists_4),
               "The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided an object that does not contain only numeric data, it includes type: character")
})

# check_pair_types
test_that("check_pair_types works", {
  expect_error(check_pair_types(data.frame(x = 1)),
               "The pair_types object")
  expect_error(check_pair_types(test_pair_types_5),
               "Your pair_types input has  ")
  expect_error(check_pair_types('test_pair_types'),
               "The pair_types object must be the output of the get_pair_types() function, but you provided:  character", fixed = TRUE)
})

#test check_subset_pairs_input works
test_that("check_subset_pairs_input works", {
  #one that works
  expect_true(is.null(check_subset_pairs_input(test_dists)))
  #one that isn't a df
  expect_error(check_subset_pairs_input("test_dists"),
               "The dists object must be a data.frame, but you provided: character",
               fixed = TRUE)
})
##################################test get_frac_intra#####################################
test_that("check_get_frac_intra_input works", {
  #normal one that works with pair_types input
  expect_null(check_get_frac_intra_input(pair_types = test_pair_types))
  #one that's input isn't pair_types input at all
  #one that's input is similar but wrong # cols
  expect_error(
    check_get_frac_intra_input(pair_types = test_pair_types_2[,1:5]),
    "The pair_types object must be the output of the get_pair_types() function, but you provided a data.frame with  5  columns.",
    fixed = TRUE
  )
  #one that's input is similar but not correct colnames
  expect_error(
    check_get_frac_intra_input(pair_types = test_pair_types_3),
    "The pair_types object must be the output of the get_pair_types() function, but the data.frame you provided has  8  columns that are not the output columns needed.",
    fixed = TRUE
  )
  #one that's input is similar but not numeric dist col
  expect_error(
    check_get_frac_intra_input(pair_types = test_pair_types_4),
    "Your pair_types input does not have numeric pairwise distances, you supplied one with type character",
    fixed = TRUE
  )
  #one with patient transfer network pair_types
  expect_null(check_get_frac_intra_input(pair_types = test_pair_types_pt_trans))
  # all null
  expect_error(check_get_frac_intra_input(pair_types = NULL),
               "Please provide the output from.")
})

##################################test get_clusters#####################################
test_that("check_get_clusters_inputs works", {
  #one that should work
  expect_true(is.null(check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = 1, bootstrap = NULL)))
  #one where tr isn't a tree
  expect_error(
    check_get_clusters_inputs(tr = "test_tr", locs = test_locs, pureness = 1, bootstrap = NULL),
    "The tr object must be a phylogenetic tree, read into R using the ape::read.tree() function. You have supplied a  character",
    fixed = TRUE
  )
  #where locs isn't a named list
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = "test_locs", pureness = 1, bootstrap = NULL),
    "The locs object must be a named list of locations named by sample IDs",
    fixed = TRUE
  )
  #where pureness isn't the right kind of value
  #not between .5 and 1
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = 100, bootstrap = NULL),
    "The pureness value must be between 0.5 and 1, you supplied a value of  100",
    fixed = TRUE
  )
  #not an int
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = "100", bootstrap = NULL),
    "The pureness value must be a double, you supplied type  character",
    fixed = TRUE
  )
  #where bootstrap isn't the right kind of value
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = 1, bootstrap = "NULL"),
    "The bootstrap value must be a double or NULL, you supplied type  character",
    fixed = TRUE
  )
  #where bootstrap isn't in the range of values
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = 1, bootstrap = 105),
    "The bootstrap value must be between 0 and 100, you supplied a value of 105",
    fixed = TRUE
  )
  #where tr is a subset of locs
  expect_warning(
    check_get_clusters_inputs(tr = test_tr_2, locs = test_locs, pureness = 1, bootstrap = NULL),
    "You have not provided the same set of locs and tip labels on the tree. Will subset",
    fixed = TRUE
  )
  #where locs is s subset of tr
  expect_warning(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs_3, pureness = 1, bootstrap = NULL),
    "You have not provided the same set of locs and tip labels on the tree. Will subset",
    fixed = TRUE
  )
  #where tr and locs have no IDs in common
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs_2, pureness = 1, bootstrap = NULL),
    "You must supply a locs object with at least two IDs in common with your tip labels of your tr object, you have provided  0  isolates in common.",
    fixed = TRUE
  )

})

# check_control_labels

test_that("check_control_labels works", {
  expect_null(check_control_labels(NULL))
  expect_error(check_control_labels(c('a','b')),
               "control_labels must be either ")
})

##################################test get_genetic_flow#####################################
test_that("check_facility_fsp_input works", {
  #one that works
  expect_true(is.null(check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_5, matrix = TRUE, pt = test_pt_5)))
  #one that works w long form
  expect_true(is.null(check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_5, matrix = FALSE, pt = test_pt_5)))
  #one that doesn't work because no more than one isolate/facility
  expect_error(
    check_facility_fsp_input(fasta = test_fasta, locs = test_locs, matrix = TRUE, pt = test_pt),
    "Please provide isolates that appear in a facility at least twice, you have not provided locs of at least two isolates in two facilities",
    fixed = TRUE
  )
  #check form input
  expect_error(
    check_facility_fsp_input(fasta = test_fasta, locs = test_locs, matrix = "hi", pt = test_pt),
    "matrix must be logical to determine output format, you have provided hi",
    fixed = TRUE
  )
  expect_error(
    check_facility_fsp_input(fasta = test_fasta, locs = test_locs, matrix = 12, pt = test_pt),
    "matrix must be logical to determine output format, you have provided 12",
    fixed = TRUE
  )
  #check fasta
  expect_error(
    check_facility_fsp_input(fasta = "test_fasta", locs = test_locs, matrix = TRUE, pt = test_pt),
    "The fasta object must be of class DNAbin, you have supplied an object of class  character",
    fixed = TRUE
  )
  #check fasta vs. locs when they don't have 2 in common
  expect_error(
    check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_6, matrix = TRUE, pt = test_pt_6),
    "Please provde a fasta object and locs object with at least two samples in common",
    fixed = TRUE
  )
  #when not two isolates in two facilities
  expect_error(
    check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_7, matrix = TRUE, pt = test_pt_7),
    "Please provide isolates that appear in a facility at least twice, you have not provided locs of at least two isolates in two facilities",
    fixed = TRUE
  )
  #check fasta vs. locs when you have to subset
  expect_warning(
    check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_9[2:6], matrix = TRUE, pt = test_pt_9[2:6]),
    "You have provided  5  isolate IDs in common between locs and your fasta. Will subset.",
    fixed = TRUE
  )
  #one where there's only isolate from a location and will subset
  expect_warning(
    check_facility_fsp_input(fasta = test_fasta_3, locs = test_locs_8, matrix = TRUE, pt = test_pt_8),
    "You have provided at least one isolate that is the only one in its location. Will subset to exclude location  1",
    fixed = TRUE
  )
})

test_that("check_long_form_input works", {
  #a good one
  expect_true(is.null(check_long_form_input(mat, c('loc1','loc2','fsp'))))
  expect_null(check_long_form_input(mat, c('a','b','c')))
  #one that is wrong type
  expect_error(is.null(check_long_form_input("mat")),
               "The fsp matrix object must be a data.frame or matrix but you provided: character",
               fixed = TRUE)
  #one not symmetric
  expect_error(is.null(check_long_form_input(mat[-1,])),
               "The fsp matrix object must be symmetric, same number of columns as rows",
               fixed = TRUE)
  #one with wrong colnames to rownames
  expect_error(is.null(check_long_form_input(mat3)),
               "The fsp matrix object must be symmetric, same number of columns as rows and columns and rows must have the same names",
               fixed = TRUE)
  expect_error(check_long_form_input(mat, c('a','b')),
               "col_names must be vector of length 3, but you provided length: 2")
  expect_error(check_long_form_input(mat, mat),
               "col_names must be a character vector but you provided: data.frame")
})

##################################test allele_freq#####################################
test_that("check_allele_freq_input works", {
  #one that works without subset
  expect_true(is.null(check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 1, alleles = c("a", "t"))))
  #one that works with subset
  expect_true(is.null(check_allele_freq_input(x = test_fasta, subset = c(TRUE, FALSE, TRUE), allele_n = 1, alleles = c("a", "t"))))
  #one with wrong x
  expect_error(
    check_allele_freq_input(x = "test_fasta", subset = NULL, allele_n = 1, alleles = c("a", "t")),
    "The x you have provided is not a DNAbin, you provided a character",
    fixed = TRUE
  )
  #one with wrong subset
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = "NULL", allele_n = 1, alleles = c("a", "t")),
    "The subset vector you have provided is not logical, you provided a character",
    fixed = TRUE
  )
  #one with wrong allele_n #
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 7, alleles = c("a", "t")),
    "The allele_n value you have provided must be 1 or 2, you have provided 7",
    fixed = TRUE
  )
  #one with wrong allele_n type
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = "1", alleles = c("a", "t")),
    "The allele_n value you have provided must be numeric 1 or 2, you have provided type character",
    fixed = TRUE
  )
  #one with one alleles
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 1, alleles = "t"),
    "The alleles vector must be a character vector of length 2, you have provided,  t",
    fixed = TRUE
  )
  #with 3 alleles
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 1, alleles = c("a","t","g")),
    "The alleles vector must be a character vector of length 2, you have provided,  a t g",
    fixed = TRUE
  )
  # expect_error(
  #   check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 1, alleles = c(1, "t")),
  #   "The alleles vector must be a character vector of length 2, you have provided,  1 t",
  #   fixed = TRUE
  # )
})

##################################test within_pop_var#####################################
test_that("check_within_pop_var_inputs works", {
  #one that works
  expect_true(is.null(check_within_pop_var_inputs(subset_snp_mat = test_fasta, subset = c(TRUE, FALSE, FALSE))))
  #one with wrong subset_snp_mat
  expect_error(
    check_within_pop_var_inputs(subset_snp_mat = test_fasta, subset = "TRUE"),
    "The subset vector you have provided is not logical, you provided a character",
    fixed = TRUE
  )
  #one with wrong subset
  expect_error(
    check_within_pop_var_inputs(subset_snp_mat = "test_fasta", subset = c(TRUE, FALSE, FALSE)),
    "The subset_snp_mat you have provided is not a DNAbin, you provided a character",
    fixed = TRUE
  )
})

##################################test find_major_alleles#####################################
test_that("check_find_major_alleles_input works", {
  #one that works (ref = null)
  expect_true(is.null(check_find_major_alleles_input(as.character(test_fasta), ref = NULL)))
  #one that works (ref != null)
  expect_true(is.null(check_find_major_alleles_input(as.character(test_fasta), ref = as.character(test_fasta)[1,])))
  #one where fasta isn't a char mat
  expect_error(
    check_find_major_alleles_input(test_fasta, ref = NULL),
    "Fasta must be a character matrix made by converting a DNAbin object using as.character, you provided DNAbin",
    fixed = TRUE
  )
  #one where ref isn't a char vec
  expect_error(
    check_find_major_alleles_input(as.character(test_fasta), ref = rep(1, ncol(test_fasta))),
    "Ref must be a character vector, you provided numeric",
    fixed = TRUE
  )
  #one where fasta and ref are different lengths
  expect_error(
    check_find_major_alleles_input(as.character(test_fasta), ref = "test"),
    "Ref must be reference sequence of the fasta file you provided, so lengths of the sequences must match. You provided a reference sequence of length 1  and a fasta of sequence length 83976",
    fixed = TRUE
  )
})


##################################test make_meta_seqs#####################################
test_that("check_make_meta_seqs_input works", {
  #one that works (pt not null)
  expect_true(is.null(check_make_meta_seqs_input(test_fasta, test_locs, test_pt)))
  #one where fasta is not DNAbin
  expect_error(
    check_make_meta_seqs_input("test_fasta", test_locs, test_pt),
    "The fasta object must be of class DNAbin, you have supplied an object of class  character",
    fixed = TRUE
  )
  #one where locs isn't a named list
  expect_error(
    check_make_meta_seqs_input(test_fasta, "test_locs", test_pt),
    "The locs object must be a named list of locations named by sample IDs",
    fixed = TRUE
  )
  #one where pt isn't a named list
  expect_error(
    check_make_meta_seqs_input(test_fasta, test_locs, "test_pt"),
    "The pt object must be a named list of patient IDs named by sample IDs or NULL",
    fixed = TRUE
  )
  #one where locs and pt don't rep the same samps
  expect_error(
    check_make_meta_seqs_input(test_fasta, test_locs[1:10], test_pt),
    "locs and pt must represent the same samples",
    fixed = TRUE
  )
  #one where locs and fasta don't rep the same samps
  expect_error(
    check_make_meta_seqs_input(test_fasta[1:10], test_locs, test_pt),
    "fasta must contain sequences for the same samples as locs and pt",
    fixed = TRUE
  )
})


##################################test reverse_list_str#####################################
test_that("check_reverse_list_str_input works", {
  #one that works
  expect_true(is.null(check_reverse_list_str_input(ls)))
  #one that isn't a list
  expect_error(
    check_reverse_list_str_input("ls"),
    "The ls object must be a list but you provided: character",
    fixed = TRUE
  )
  #one that is a list but not of lists
  expect_error(
    check_reverse_list_str_input(ls_2),
    "The ls object must be a list of lists but you provided at lease one element that is not a list",
    fixed = TRUE
  )
})

##################################test get_largest_subtree#####################################
test_that("check_get_largest_subtree_input works", {
  #one that works
  expect_true(is.null(check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = 1)))
  #one where subtrees is wrong
  #not right
  expect_error(
    check_get_largest_subtree_input(subtrs = "test_subtr", isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The subtrs object must be the result of calling ape::subtrees, the object you have provided is not a list",
    fixed = TRUE
  )
  #list but not phylo
  expect_error(
    check_get_largest_subtree_input(subtrs = ls, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The subtrs object must be the result of calling ape::subtrees, the contents of the list you provided are not of type phylo",
    fixed = TRUE
  )
  #one where isolate labels is wrong
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = "test_pt", control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The locs object must be a named list of locations named by sample IDs",
    fixed = TRUE
  )
  #one where they don't match at all
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt_2, control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The subtrs object must include all of the isolates provided in the isolate_labels object",
    fixed = TRUE
  )
  #one where locs is a subset of tr
  expect_true(is.null(check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt_3, control_labels = NULL, bootstrap = NULL, pureness = 1)))
  #one where tree is a subset of locs
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr_2, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The subtrs object must include all of the isolates provided in the isolate_labels object",
    fixed = TRUE
  )
  #one where control labels is wrong
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = c('A'='A'), bootstrap = NULL, pureness = 1),
    "The subtrs object must include all of the isolates provided in the control_labels object",
    fixed = TRUE
  )
  #one where pureness is wrong
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = "1"),
    "The pureness value must be a double, you supplied type  character",
    fixed = TRUE
  )
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = 100),
    "The pureness value must be between 0.5 and 1, you supplied a value of  100",
    fixed = TRUE
  )
  #one where bootstrap is wrong
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = "NULL", pureness = 1),
    "The bootstrap value must be a double or NULL, you supplied type  character",
    fixed = TRUE
  )
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = 1000, pureness = 1),
    "The bootstrap value must be between 0 and 100, you supplied a value of 1000",
    fixed = TRUE
  )
})



##################################test patient_flow#####################################
test_that("check_get_patient_flow_input works", {
  #one that works with pair_types
  expect_null(check_get_patient_flow_input(pt_trans_df = test_pt_trans_net, paths = FALSE, locs = NULL))
  #one where pt_trans_df is wrong
  expect_error(check_get_patient_flow_input(pt_trans_df = "test_pt_trans_net", paths = FALSE, locs = NULL),
               "The pt_trans_df object must be a data.frame or matrix, you provided a  character",
               fixed = TRUE)
  #one where pt_trans_df has wrong ncol
  expect_error(check_get_patient_flow_input(pt_trans_df = test_pt_trans_net_2, paths = FALSE, locs = NULL),
               "The pt_trans_df object must be a data.frame or matrix with 3 columns, you provided  2",
               fixed = TRUE)
  #one where pt_trans_df has wrong colnames
  expect_error(check_get_patient_flow_input(pt_trans_df = test_pt_trans_net_3, paths = FALSE, locs = NULL),
               "The pt_trans_df object must be a data.frame or matrix with 3 columns named 'source_facil', 'dest_facil', and 'n_transfers', you provided  A B C",
               fixed = TRUE)
  #one where pt_trans_df has wrong coltypes
  expect_error(check_get_patient_flow_input(pt_trans_df = test_pt_trans_net_4, paths = FALSE, locs = NULL),
               "The pt_trans_df object must be a data.frame or matrix with 3 columns named 'source_facil', 'dest_facil', and 'n_transfers', of types character, character and numeric consecutively. you provided  factor factor character",
               fixed = TRUE)
  #one where paths is wrong
  expect_error(check_get_patient_flow_input(pt_trans_df = test_pt_trans_net_5, paths = 12, locs = NULL),
               "paths argument must be a logical/logical value representing whether you want to return the shortest paths used to generate the indirect flow metric. You have provided numeric",
               fixed = TRUE)

  expect_error(check_get_patient_flow_input(pt_trans_df = dplyr::bind_rows(test_pt_trans_net, test_pt_trans_net), paths = FALSE, locs = NULL),
               "Multiple rows in the patient transfer network contain the same source and destination facility. Please include only unique source and destination pairs.",
               fixed = TRUE)
})

# check_thresh

test_that("check_thresh works", {
  expect_null(check_thresh(1))
  expect_error(check_thresh(c(1,2)),
               "thresh must be a single value, you provided")
  expect_error(check_thresh('a'),
               "thresh must be a numeric value, you provided character")
  expect_error(check_thresh(-1),
               "thresh must be a positive numeric value, you provided")
})

# check_paths

test_that("check_paths works", {
  expect_null(check_paths(TRUE))
  expect_error(check_paths(1),
               "paths argument must be a logical/logical value representing whether you want to return the shortest paths used to generate the indirect flow metric. You have provided numeric")
})


# check_summarize_pairs_input

test_that("check_summarize_pairs_input works", {
  expect_null(check_summarize_pairs_input(test_pair_types, 'min', 5))
  expect_null(check_summarize_pairs_input(test_pair_types, c('min', 'median'), 5))
  expect_null(check_summarize_pairs_input(test_pair_types, 'min', c(5, 10)))
  expect_null(check_summarize_pairs_input(test_pair_types, NULL, 5))
  expect_null(check_summarize_pairs_input(test_pair_types, 'min', NULL))
  expect_error(check_summarize_pairs_input('test_pair_types', min, 5),
               "The pair_types object must be the output of the get_pair_types() function, but you provided:  character", fixed = TRUE)
  expect_error(check_summarize_pairs_input(test_pair_types, min, 5),
               "The summary_fns argment must either be `NULL` or a character vector of function names you wish to use to summarize inter-facility pariwise distances. You have provided function")
  expect_error(check_summarize_pairs_input(test_pair_types, c('min','asd'), 5),
               "The summary_fns argment must either be `NULL` or a character vector of function names you wish to use to summarize inter-facility pairwise distances. You have provided at least one element that is not a function:")
  expect_error(check_summarize_pairs_input(test_pair_types, 'min', '5'),
               "The threshs argment must either be `NULL` or a numeric vector of the pairwise SNV distance thresholds you wish to use to summarize inter-facility pairs. You have provided character")

})

# check_merge_inter_summaries_input

test_that("check_merge_inter_summaries_input works", {
  expect_null(check_merge_inter_summaries_input(test_pt_flow, test_isolate_pair_summary, test_fsp_long))
  expect_null(check_merge_inter_summaries_input(test_pt_flow, test_isolate_pair_summary, NULL))
  expect_null(check_merge_inter_summaries_input(test_pt_flow, NULL, test_fsp_long))
  expect_null(check_merge_inter_summaries_input(NULL, test_isolate_pair_summary, test_fsp_long))
  expect_null(check_merge_inter_summaries_input(NULL, NULL, NULL))
  expect_error(check_merge_inter_summaries_input('test_pt_flow', test_isolate_pair_summary, test_fsp_long),
               "patient_flow must be a data.frame but you provided character")
  expect_error(check_merge_inter_summaries_input(test_pt_flow, 'test_isolate_pair_summary', test_fsp_long),
               "isolate_pair_summary must be a data.frame but you provided character")
  expect_error(check_merge_inter_summaries_input(test_pt_flow, test_isolate_pair_summary, 'test_fsp_long'),
               "fsp_long must be a data.frame but you provided character")
  expect_error(check_merge_inter_summaries_input(test_pt_flow[,3:5], test_isolate_pair_summary, test_fsp_long),
               "There must be at least one column in common between all inputs.")
})
