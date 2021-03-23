#tests for regentrans - to be made when the public facing dataset is available
####################################test data######################################
test_locs <- locs[1:4]
test_locs_2 <- locs[5:8]
test_locs_3 <- locs[1:3]
test_locs_4 <- locs[1]
test_locs_5 <- locs[locs %in% c("A", "F", "H")]
test_locs_6 <- locs[1:4]
test_locs_7 <- test_locs_5[1:4]
test_locs_8 <- test_locs_5[2:6]
test_locs_9 <- test_locs_5[1:6]
test_pt <- pt[1:4]
test_pt_2 <- pt[5:8]
test_pt_3 <- pt[1:3]
test_dists <- dists[names(test_locs), names(test_locs)]
test_dists_2 <- dists[names(test_locs_3), names(test_locs_3)]
test_snv_dists <- get_snv_dists(dists = test_dists, locs = test_locs, pt = test_pt)
test_snv_dists_2 <- test_snv_dists[,2:ncol(test_snv_dists)]
test_snv_dists_3 <- test_snv_dists
colnames(test_snv_dists_3) <- c("A", "B", "C", "D", "E", "F", "G", "H")
test_snv_dists_4 <- test_snv_dists
test_snv_dists_4$Pairwise_Dists <- as.character(test_snv_dists_4$Pairwise_Dists)
test_tr <- keep.tip(tr,names(test_pt))
test_tr_2 <- keep.tip(tr,names(test_pt_3))
test_fasta <- fasta[names(test_locs),]
test_fasta_2 <- fasta[names(test_locs_5),]
test_fasta_3 <- fasta[names(test_locs_9),]
ls <- list(list(a = 2, b = 3), list(c = "a", d = "b"))
ls_2 <- list(list(a = 2, b = 3), list(c = "a", d = "b"), "x")
test_subtr <- ape::subtrees(test_tr)
test_subtr_2 <- ape::subtrees(test_tr_2)

mat <- data.frame(matrix(data = c(0, 20, 12,
                                  20, 0, 26,
                                  30, 26, 0), nrow = 3, ncol = 3))
rownames(mat) <- c("A", "B", "C")
colnames(mat) <- c("A", "B", "C")
test_pt_trans_net <- na.omit(data.frame(as.table(as.matrix(mat))))
#pat_flow <- dplyr::bind_cols(pat_flow %>% filter(Var1 != Var2))
colnames(test_pt_trans_net) <- c("source_facil", "dest_facil", "n_transfers")
test_pt_trans_net$n_transfers <- as.numeric(test_pt_trans_net$n_transfers)
test_pt_trans_net_2 <- test_pt_trans_net[,2:ncol(test_pt_trans_net)]
test_pt_trans_net_3 <- test_pt_trans_net
colnames(test_pt_trans_net_3) <- c("A", "B", "C")
test_pt_trans_net_4 <- test_pt_trans_net
test_pt_trans_net_4$n_transfers <- as.character(test_pt_trans_net_4$n_transfers)


####################################test get_snv_dists######################################
test_that("check_get_snv_dists_input works", {
  #check ones that should pass
  #with pt
  expect_true(is.null(check_get_snv_dists_input(test_dists, test_locs, test_pt)))
  #without pt
  expect_true(is.null(check_get_snv_dists_input(test_dists, test_locs, pt=NULL)))
  #one input - dists
  expect_error(
    check_get_snv_dists_input(dists = test_dists),
    'argument "locs" is missing, with no default'
  )
  #one input - locs
  expect_error(
    check_get_snv_dists_input(locs = test_locs),
    'argument "dists" is missing, with no default'
  )
  #no locs input
  expect_error(
    check_get_snv_dists_input(dists = test_dists, pt = test_pt),
    'argument "locs" is missing, with no default'
  )
  #test_dists not a dists object
  expect_error(
    check_get_snv_dists_input(dists = "test_dists", locs = test_locs),
    "The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided: character"
  )
  #locs not a named list
  expect_error(
    check_get_snv_dists_input(dists = test_dists, locs = "test_locs"),
    "The locs object must be a a named list of locations named by sample IDs"
  )
  #pt not a named list
  expect_error(
    check_get_snv_dists_input(dists = test_dists, locs = test_locs, pt = "test_pt"),
    "The pt object must be a a named list of locations named by sample IDs"
  )
  #dist names not matching locs names
  expect_error(
    check_get_snv_dists_input(dists = test_dists, locs = test_locs_2, pt = test_pt),
    "You have not provided locations of at least 2 isolates in your SNV distance matrix (dists). Please provide locations for at least 2 isolates in your SNV distance matrix."
  )
  #dist names not matching pt names
  expect_error(
    check_get_snv_dists_input(dists = test_dists, locs = test_locs, pt = test_pt_2),
    "You have not provided patient IDs for at least 2 isolates in your SNV distance matrix (dists). Please provide patient IDs for at least 2 isolates in your SNV distance matrix."
  )
  #one with less than two in common
  expect_error(
    check_get_snv_dists_input(dists = test_dists, locs = test_locs_4, pt = test_pt),
    "You have only supplied locations for 1 isolates. Please supply a named vector of locations for at least 2 isolates"
  )
  #warnings
  #dist names being a subset of locs names
  expect_warning(
    check_get_snv_dists_input(dists = test_dists_2, locs = test_locs, pt = test_pt),
    "1: In check_dists_vs_locs(dists, locs) :
  You have supplied a list of more isolates (n =  4 ) with locations than exist in your SNV distance matrix (n =  3 ). Will subset
2: In check_dists_vs_locs(dists, locs) :
  You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset"
  )
  #locs names being a subset of dists names
  expect_warning(
    check_get_snv_dists_input(dists = test_dists, locs = test_locs_3, pt = test_pt),
    "1: In check_dists_vs_locs(dists, locs) :
  You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset
2: In check_pt_vs_locs(pt, locs) :
  You have supplied a patient vector of length  4  and  location vector of length  3 . We will subset these lists so that they have the same isolates.
3: In check_pt_vs_locs(pt, locs) :
  You have not supplied patient IDs (pt) and locations (locs) for the same samples. We will these vectors so that they contain the same isolates."
  )
  #pt names being a subset of dists names
  expect_warning(
    check_get_snv_dists_input(dists = test_dists, locs = test_locs, pt = test_pt_3),
    "1: In check_pt_vs_locs(pt, locs) :
  You have supplied a patient vector of length  3  and  location vector of length  4 . We will subset these lists so that they have the same isolates.
2: In check_pt_vs_locs(pt, locs) :
  You have not supplied patient IDs (pt) and locations (locs) for the same samples. We will these vectors so that they contain the same isolates."
  )
  #dists names being a subset of pt names
  expect_warning(
    check_get_snv_dists_input(dists = test_dists_2, locs = test_locs, pt = test_pt),
    "1: In check_dists_vs_locs(dists, locs) :
  You have supplied a list of more isolates (n =  4 ) with locations than exist in your SNV distance matrix (n =  3 ). Will subset
2: In check_dists_vs_locs(dists, locs) :
  You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset"
  )
})

##################################test get_frac_intra#####################################
test_that("check_get_frac_intra_input works", {
  #normal one that works with snv_dists input
  expect_false(check_get_frac_intra_input(snv_dists = test_snv_dists, threshs = seq(1,50,1), dists = NULL, locs = NULL, pt = NULL))
  #normal one that works without snv_dists input
  expect_true(check_get_frac_intra_input(snv_dists = NULL, threshs = seq(1,50,1), dists = test_dists, locs = test_locs, pt = test_pt))
  #one that's input isn't snv_dists input at all
  expect_error(
    check_get_frac_intra_input(snv_dists = "test_snv_dists", threshs = seq(1,50,1), dists = NULL, locs = NULL, pt = NULL),
    "The snv_dists object must be the output of the get_snv_dists() function, but you provided:  character"
  )
  #one that's input is similar but wrong # cols
  expect_error(
    check_get_frac_intra_input(snv_dists = test_snv_dists_2, threshs = seq(1,50,1), dists = NULL, locs = NULL, pt = NULL),
    "The snv_dists object must be the output of the get_snv_dists() function, but you provided a data.frame with  7  columns."
  )
  #one that's input is similar but not correct rownames
  expect_error(
    check_get_frac_intra_input(snv_dists = test_snv_dists_3, threshs = seq(1,50,1), dists = NULL, locs = NULL, pt = NULL),
    "The snv_dists object must be the output of the get_snv_dists() function, but the data.frame you provided has  8  columns that are not the output columns needed."
  )
  #one that's input is similar but not numeric dist col
  expect_error(
    check_get_frac_intra_input(snv_dists = test_snv_dists_4, threshs = seq(1,50,1), dists = NULL, locs = NULL, pt = NULL),
    "Your snv_dists input does not have numeric pairwise distances, you supplied one with type character"
  )
  #one that's input doesn't work for get_snv_dists
  expect_error(
    check_get_frac_intra_input(snv_dists = NULL, threshs = seq(1,50,1), dists = "test_dists", locs = test_locs, pt = test_pt),
    "The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided: character"
  )
  #one where threshs isn't a numeric vector
  expect_error(
    check_get_frac_intra_input(snv_dists = test_snv_dists, threshs = "seq(1,50,1)", dists = NULL, locs = NULL, pt = NULL),
    "threshs must be a numeric vector, you provided a  character"
  )
})

##################################test get_clusters#####################################
test_that("check_get_clusters_inputs works", {
  #one that should work
  expect_true(is.null(check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = 1, bootstrap = NULL)))
  #one where tr isn't a tree
  expect_error(
    check_get_clusters_inputs(tr = "test_tr", locs = test_locs, pureness = 1, bootstrap = NULL),
    "The tr object must be a phylogenetic tree, read into R using the ape::read.tree() function. You have supplied a  character"
  )
  #where locs isn't a named list
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = "test_locs", pureness = 1, bootstrap = NULL),
    "The locs object must be a a named list of locations named by sample IDs"
  )
  #where pureness isn't the right kind of value
  #not between .5 and 1
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = 100, bootstrap = NULL),
    "The pureness value must be between 0.5 and 1, you supplied a value of  100"
  )
  #not an int
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = "100", bootstrap = NULL),
    "The pureness value must be a double, you supplied type  character"
  )
  #where bootstrap isn't the right kind of value
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = 1, bootstrap = "NULL"),
    "The bootstrap value must be a double or NULL, you supplied type  character"
  )
  #where bootstrap isn't in the range of values
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs, pureness = 1, bootstrap = 105),
    "The bootstrap value must be between 0 and 100, you supplied a value of 105"
  )
  #where tr is a subset of locs
  expect_warning(
    check_get_clusters_inputs(tr = test_tr_2, locs = test_locs, pureness = 1, bootstrap = NULL),
    "You have not provided the same set of locs and tip labels on the tree. Will subset"
  )
  #where locs is s subset of tr
  expect_warning(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs_3, pureness = 1, bootstrap = NULL),
    "You have not provided the same set of locs and tip labels on the tree. Will subset"
  )
  #where tr and locs have no IDs in common
  expect_error(
    check_get_clusters_inputs(tr = test_tr, locs = test_locs_2, pureness = 1, bootstrap = NULL),
    "You must supply a locs object with at least two IDs in common with your tip labels of your tr object, you have provided  0  isolates in common."
  )

})

##################################test get_facility_fsp#####################################
test_that("check_facility_fsp_input works", {
  #one that works
  expect_true(is.null(check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_5, form = "matrix")))
  #one that works w long form
  expect_true(is.null(check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_5, form = "long")))
  #one that doesn't work because no more than one isolate/facility
  expect_error(
    check_facility_fsp_input(fasta = test_fasta, locs = test_locs, form = "matrix"),
    "Please provide isolates that appear hospitals at least once"
  )
  #check form input
  expect_error(
    check_facility_fsp_input(fasta = test_fasta, locs = test_locs, form = "hi"),
    "form must be either 'long' or 'matrix' to determine output format, you have provided hi"
  )
  expect_error(
    check_facility_fsp_input(fasta = test_fasta, locs = test_locs, form = 12),
    "form must be either 'long' or 'matrix' to determine output format, you have provided 12"
  )
  #check fasta
  expect_error(
    check_facility_fsp_input(fasta = "test_fasta", locs = test_locs, form = "matrix"),
    "The fasta object must be of class DNAbin, you have supplied an object of class  character"
  )
  #check fasta vs. locs when they don't have 2 in common
  expect_error(
    check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_6, form = "matrix"),
    "Please provde a fasta object and locs object with at least two samples in common"
  )
  #when not two isolates in two facilities
  expect_error(
    check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_7, form = "matrix"),
    "Please provide isolates that appear in a facility at least twice, you have not provided locs of at least two isolates in two facilities"
  )
  #check fasta vs. locs when you have to subset
  expect_warning(
    check_facility_fsp_input(fasta = test_fasta_2, locs = test_locs_8, form = "matrix"),
    "You have provided  5  isolate IDs in common between locs and your fasta. Will subset."
  )
  #one where there's only isolate from a location and will subset
  expect_warning(
    check_facility_fsp_input(fasta = test_fasta_3, locs = test_locs_9, form = "matrix"),
    "You have provided at least one isolate that is the only one in its location. Will subset to exclude location  1"
  )
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
    "The x you have provided is not a DNAbin, you provided a character"
  )
  #one with wrong subset
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = "NULL", allele_n = 1, alleles = c("a", "t")),
    "The subset vector you have provided is not logical, you provided a character"
  )
  #one with wrong allele_n #
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 7, alleles = c("a", "t")),
    "The allele_n value you have provided must be 1 or 2, you have provided 7"
  )
  #one with wrong allele_n type
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = "1", alleles = c("a", "t")),
    "The allele_n value you have provided must be numeric 1 or 2, you have provided type character"
  )
  #one with one alleles
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 1, alleles = "t"),
    "The alleles vector must be a character vector of length 2, you have provided,  t"
  )
  #with 3 alleles
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 1, alleles = c("a","t","g")),
    "The alleles vector must be a character vector of length 2, you have provided,  a t g"
  )
  expect_error(
    check_allele_freq_input(x = test_fasta, subset = NULL, allele_n = 1, alleles = c(1, "t")),
    "The alleles vector must be a character vector of length 2, you have provided,  1 t"
  )
})

##################################test within_pop_var#####################################
test_that("within_pop_var_input_checks works", {
  #one that works
  expect_true(is.null(within_pop_var_input_checks(subset_snp_mat = test_fasta, subset = c(TRUE, FALSE, FALSE))))
  #one with wrong subset_snp_mat
  expect_error(
    within_pop_var_input_checks(subset_snp_mat = test_fasta, subset = "TRUE"),
    "The subset vector you have provided is not logical, you provided a character"
  )
  #one with wrong subset
  expect_error(
    within_pop_var_input_checks(subset_snp_mat = "test_fasta", subset = c(TRUE, FALSE, FALSE)),
    "The subset_snp_mat you have provided is not a DNAbin, you provided acharacter"
  )
})

##################################test reverse_list_str#####################################
test_that("check_reverse_list_str_input works", {
  #one that works
  expect_true(is.null(check_reverse_list_str_input(ls)))
  #one that isn't a list
  expect_error(
    check_reverse_list_str_input("ls"),
    "The ls object must be a list but you provided: character"
  )
  #one that is a list but not of lists
  expect_error(
    check_reverse_list_str_input(ls_2),
    "The ls object must be a list of lists but you provided at lease one element that is not a list"
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
    "The subtrs object must be the result of calling ape::subtrees, the object you have provided is not a list"
  )
  #list but not phylo
  expect_error(
    check_get_largest_subtree_input(subtrs = ls, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The subtrs object must be the result of calling ape::subtrees, the contents of the list you provided are not of type phylo"
  )
  #one where isolate labels is wrong
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = "test_pt", control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The locs object must be a a named list of locations named by sample IDs"
  )
  #one where they don't match at all
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt_2, control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The subtrs object must include all of the isolates provided in the isolate_labels object"
  )
  #one where locs is a subset of tr
  expect_true(is.null(check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt_3, control_labels = NULL, bootstrap = NULL, pureness = 1)))
  #one where tree is a subset of locs
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr_2, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = 1),
    "The subtrs object must include all of the isolates provided in the,  isolate_labels object"
  )
  #one where control labels is wrong
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = "NULL", bootstrap = NULL, pureness = 1),
    "The subtrs object must include all of the isolates provided in the,  control_labels object"
  )
  #one where pureness is wrong
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = "1"),
    "The pureness value must be a double, you supplied type  character"
  )
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = NULL, pureness = 100),
    "The pureness value must be between 0.5 and 1, you supplied a value of  100"
  )
  #one where bootstrap is wrong
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = "NULL", pureness = 1),
    "The bootstrap value must be a double or NULL, you supplied type  character"
  )
  expect_error(
    check_get_largest_subtree_input(subtrs = test_subtr, isolate_labels = test_pt, control_labels = NULL, bootstrap = 1000, pureness = 1),
    "The bootstrap value must be between 0 and 100, you supplied a value of 1000"
  )
})




##################################test patient_flow#####################################
test_that("check_pt_transfer_input works", {
  #one that works with snv_dists
  expect_false(check_pt_transfer_input(pt_trans_net = test_pt_trans_net, snv_dists = test_snv_dists, thresh = 50, dists = NULL, locs = NULL, pt = NULL))
  #one where we make snv_dists
  expect_true(check_pt_transfer_input(pt_trans_net = test_pt_trans_net, snv_dists = NULL, thresh = 50, dists = test_dists, locs = test_locs, pt = test_pt))
  #one with wrong snv_dists input
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net, snv_dists = "test_snv_dists", thresh = 50, dists = NULL, locs = NULL, pt = NULL),
               "The snv_dists object must be the output of the get_snv_dists() function, but you provided:  character")
  #one with similar to snv_dists but not right ncols
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net, snv_dists = test_snv_dists_2, thresh = 50, dists = NULL, locs = NULL, pt = NULL),
               "The snv_dists object must be the output of the get_snv_dists() function, but you provided a data.frame with  7  columns.")
  #one similar but not right colnames
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net, snv_dists = test_snv_dists_3, thresh = 50, dists = NULL, locs = NULL, pt = NULL),
               "The snv_dists object must be the output of the get_snv_dists() function, but the data.frame you provided has  8  columns that are not the output columns needed.")
  #one similar with wrong coltype
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net, snv_dists = test_snv_dists_4, thresh = 50, dists = NULL, locs = NULL, pt = NULL),
               "Your snv_dists input does not have numeric pairwise distances, you supplied one with type character")
  #one with wrong input to get_snv_dists
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net, snv_dists = NULL, thresh = 50, dists = "test_dists", locs = test_locs, pt = test_pt),
              "The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided: character")
  #one where thresh is wrong
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net, snv_dists = test_snv_dists, thresh = -1, dists = NULL, locs = NULL, pt = NULL),
               "thresh must be a positive numeric value, you provided -1")
  #one where pt_trans_net is wrong
  expect_error(check_pt_transfer_input(pt_trans_net = "test_pt_trans_net", snv_dists = test_snv_dists, thresh = 50, dists = NULL, locs = NULL, pt = NULL),
               "The pt_trans_net object must be a data.frame or matrix, you provided a  character")
  #one where pt_trans_net has wrong ncol
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net_2, snv_dists = test_snv_dists, thresh = 50, dists = NULL, locs = NULL, pt = NULL),
               "The pt_trans_net object must be a data.frame or matrix with 3 columns, you provided  2")
  #one where pt_trans_net has wrong colnames
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net_3, snv_dists = test_snv_dists, thresh = 50, dists = NULL, locs = NULL, pt = NULL),
               "The pt_trans_net object must be a data.frame or matrix with 3 columns named 'source_facil', 'dest_facil', and 'n_transfers', you provided  A B C")
  #one where pt_trans_net has wrong coltypes
  expect_error(check_pt_transfer_input(pt_trans_net = test_pt_trans_net_4, snv_dists = test_snv_dists, thresh = 50, dists = NULL, locs = NULL, pt = NULL),
               "The pt_trans_net object must be a data.frame or matrix with 3 columns named 'source_facil', 'dest_facil', and 'n_transfers', of types character, character and numeric consecutively. you provided  factor factor character ")


})
