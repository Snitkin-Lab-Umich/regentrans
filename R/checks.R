#checks for regentrans inputs
#########################################################get_snv_dists####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR get_snv_dists FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
#checks dists input to get_snv_dists function
check_dists <- function(dists){
  #if it is not a matrix
  if(!(any(class(dists) == "matrix") || any(class(dists) == "data.frame"))){
    stop(paste("The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided:",
               class(dists)))
  }
  #make sure it is a correlation matrix (aka equal dimensions)
  if(dim(dists)[1] != dim(dists)[2]){
    stop(paste("The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but the dimensions of your matrix are not equal. The matrix you provided is",
               dim(dists)[1], "x", dim(dists)[2]))
  }
  #make sure there are at least 2 columns and 2 rows? Don't just want to correlate one sample to itself
  #only need ot check one dimension because the previous check made sure they were the same
  if(nrow(dists) < 2){
    stop(paste("Your SNV matrix only has ", nrow(dists), " samples. Please use a SNV distance matrix that includes 2 or more samples"))
  }
  #check that the data is numeric
  if(!(all(sapply(dists, class) == "numeric") || all(sapply(dists, class) == "integer"))){
    stop(paste("The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided an object that does not contain only numeric data, it includes type:",
               unique(sapply(dists, class))))
  }
  #check that all of the data is >= 0
  if(any(dists < 0)){
    stop("The dists object must be a SNV distance matrix returned by the dist.dna function from the ape package, but you provided an object with values < 0")
  }

}

#checks locs input to get_snv_dists function
check_locs <- function(locs){
  #check that the locs object is a named vector
  #if it is not a vector or has no names then not good (check both of these in tests)
  if(!(is.vector(locs)) || is.null(names(locs))){
    stop("The locs object must be a a named list of locations named by sample IDs")
  }
  #check that the locs object has at least 2 items
  if(length(locs) < 2){
    stop(paste("You have only supplied locations for "), length(locs),
         " isolates. Please supply a named vector of locations for at least 2 isolates")
  }
}

#checks pt input to get_snv_dists function
check_pt <- function(pt, dists){
  #check that the pt object is a named vector
  #if it is not a vector or has no names then not good (check both of these in tests)
  if(!(is.vector(pt)) | is.null(names(pt))){
    stop("The pt object must be a a named list of locations named by sample IDs")
  }
  #check that the pt object has at least 2 items
  if(length(pt) < 2){
    stop(paste("You have only supplied patient IDs for "), length(pt),
         " isolate. Please supply a named vector of patient IDs for at least 2 isolates")
  }
  #check that there are at least 2 dists and locs in common
  if(length(intersect(rownames(dists), names(pt))) < 2){
    stop(paste("You have not provided patient IDs for at least 2 isolates in your SNV distance matrix (dists). Please provide patient IDs for at least 2 isolates in your SNV distance matrix."))
  }
}

#check that samples and lengths of pt and locs vectors are the same
check_pt_vs_locs <- function(pt, locs){
  #check that there are the same number of samples in the pt vector as in the locs vector
  #where would I want to subet myself?
  if(length(pt) != length(locs)){
    warning(paste("You have supplied a patient vector of length ", length(pt), " and  location vector of length ", length(locs), ". We will subset these lists so that they have the same isolates."))
  }
  #check that the names match between pt and locs
  if(!setequal(names(pt), names(locs))){
    warning("You have not supplied patient IDs (pt) and locations (locs) for the same samples. We will subset these vectors so that they contain the same isolates.")
  }
}

#check that the names of the isolates in locs actually exist in the SNV matrix
check_dists_vs_locs <- function(dists, locs){
  #check that there are less than or equal to the number of samples in the vector than in the dists matrix
  if(length(locs) > nrow(dists)){
    warning(paste("You have supplied a list of more isolates (n = ", length(locs),
               ") with locations than exist in your SNV distance matrix (n = ",
               nrow(dists),
               "). Will subset"))
  }
  #check that there are at least 2 dists and locs in common
  if(length(intersect(rownames(dists), names(locs))) < 2){
    stop(paste("You have not provided locations of at least 2 isolates in your SNV distance matrix (dists). Please provide locations for at least 2 isolates in your SNV distance matrix."))
  }
  #warn if they will be subsetting??
  if(!setequal(names(locs), rownames(dists))){
    warning("You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset")
  }
}

check_pt_trans_net <- function(pt_trans_net, locs){
  #default is null
  if(!is.null(pt_trans_net)){
    #make sure it is a dataframe of three columns, source, dest, n_transfers
    if(!(any(class(pt_trans_net) == "data.frame") || any(class(pt_trans_net) == "matrix"))){
      stop(paste("The pt_trans_net object must be a data.frame or matrix, you provided a ",
                 class(pt_trans_net)))
    }
    #make sure three cols
    if(ncol(pt_trans_net) != 3){
      stop(paste("The pt_trans_net object must be a data.frame or matrix with 3 columns, you provided ",
                 ncol(pt_trans_net)))
    }
    #make sure cols have correct names
    if(!all(colnames(pt_trans_net) == c("source_facil", "dest_facil", "n_transfers" ))){
      stop(paste("The pt_trans_net object must be a data.frame or matrix with 3 columns named 'source_facil', 'dest_facil', and 'n_transfers', you provided ",
                 paste(colnames(pt_trans_net), sep = " ", collapse = " ")))
    }
    #make sure the column types are good
    if(!(all(lapply(pt_trans_net, class) == c("factor", "factor", "numeric")) || all(lapply(pt_trans_net, class) == c("character", "character", "numeric")))){
      stop(paste("The pt_trans_net object must be a data.frame or matrix with 3 columns named 'source_facil', 'dest_facil', and 'n_transfers', of types character, character and numeric consecutively. you provided ",
                 paste(lapply(pt_trans_net, class), sep = " ", collapse = " ")))
    }
    #make sure facilities have something in common with the locs and say we will subset if they are not all in common
    #make sure there are at least two in common
    if(as.numeric(length(intersect(unique(c(as.character(pt_trans_net$source_facil), as.character(pt_trans_net$dest_facil))), unique(locs)))) < 2){
      stop(paste("The pt_trans_net don't have at least two locations in common with the locs object. "))
    }
    #if not all in common warn you will subset
    !setequal(unique(c(as.character(pt_trans_net$source_facil), as.character(pt_trans_net$dest_facil))), unique(locs))
    if(!setequal(unique(c(as.character(pt_trans_net$source_facil), as.character(pt_trans_net$dest_facil))), unique(locs))){
      warning(paste("Not all of the locations you have provided between locs and the pt_trans_network match. Will subset. "))
    }
    if(length(paste(pt_trans_net$source_facil, pt_trans_net$dest_facil)) != length(unique(paste(pt_trans_net$source_facil, pt_trans_net$dest_facil)))){
      stop(paste("Multiple rows in the patient transfer network contain the same source and destination facility. Please include only unique source and destination pairs."))
    }
  }
}


#wrapper for input to snv_dists
check_get_snv_dists_input <- function(dists, locs, pt, pt_trans_net){
  #check everything that is common (aka no pt) first
  #check that the dists object is the snv object returned by dist.dna
  check_dists(dists)
  #check that the locs object is a named vector
  check_locs(locs)
  #check that the locs names exist in the dists dataframe
  check_dists_vs_locs(dists, locs)
  #if there is a pt
  if(!is.null(pt)){
    #check that the pt object is a named vector
    check_pt(pt, dists)
    #check that the pt and dists objects have the same lengths
    check_pt_vs_locs(pt, locs)
  }
  #check whether pt_trans_net is null or a directed graph
  if(!is.null(pt_trans_net)){
    check_pt_trans_net(pt_trans_net, locs)
  }
}

#check subset pairs input
check_subset_pairs_input <- function(snv_dists){
  if(!(any(class(snv_dists) == "matrix") || any(class(snv_dists) == "data.frame"))){
    stop(paste("The dists object must be a data.frame, but you provided:",
               class(snv_dists)))
  }
}

#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_snv_dists FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
#########################################################get_frac_intra####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR get_frac_intra FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
check_snv_dists <- function(snv_dists){
  #check the type
  if(!(class(snv_dists) == "data.frame")){
    stop(paste("The snv_dists object must be the output of the get_snv_dists() function, but you provided: ",
               class(snv_dists)))
  }
  #check the number of columns
  if(!(ncol(snv_dists) %in% 6:12)){
    stop(paste("The snv_dists object must be the output of the get_snv_dists() function, but you provided a data.frame with ",
               ncol(snv_dists), " columns."))
  }
  #check the colnames
  if(!((ncol(snv_dists) == 8 &&
        all(colnames(snv_dists) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2", "Patient1",  "Patient2", "Pair_Type"))) ||
       (ncol(snv_dists) == 12 &&
        all(colnames(snv_dists) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2", "Patient1",  "Patient2", "Pair_Type", "n_1_to_2_transfers", "n_2_to_1_transfers", "indirect_flow_metric_1_to_2", "indirect_flow_metric_2_to_1"))) ||
       (ncol(snv_dists) == 10 &&
        all(colnames(snv_dists) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2", "Pair_Type", "n_1_to_2_transfers", "n_2_to_1_transfers", "indirect_flow_metric_1_to_2", "indirect_flow_metric_2_to_1"))) ||
       (ncol(snv_dists) == 6 &&
        all(colnames(snv_dists) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2", "Pair_Type"))))){
    stop(paste("The snv_dists object must be the output of the get_snv_dists() function, but the data.frame you provided has ",
               ncol(snv_dists), " columns that are not the output columns needed."))
  }
  #check that there is at least one row
  if(nrow(snv_dists) < 1){
      stop(paste("Your snv_dists input has ", nrow(snv_dists), " facility pairs. Please use an snv_dists input that has 1 or more pairs (rows)"))
  }
  #check pairwise dist column is numeric
  if(!class(snv_dists$Pairwise_Dists) == "numeric"){
    stop(paste("Your snv_dists input does not have numeric pairwise distances, you supplied one with type",
               class(snv_dists$Pairwise_Dists)))
  }

}

check_get_frac_intra_input <- function(snv_dists, dists, locs, pt, pt_trans_net){
  #if SNV_dists doesnt exist and they didn't input locs and dists
  if(is.null(snv_dists) & is.null(dists) & is.null(locs)){
    stop("Please provide either an SNV dists matrix or dists and locs objecst so we can generate one for you.")
  }
  #make the SNV dists object if it needs to be made
  if(is.null(snv_dists)){
    #checks for get_snv_dists are done within here
    snv_dists <- get_snv_dists(dists = dists, locs = locs, pt = pt, pt_trans_net = pt_trans_net)
    run_snv_dists <- TRUE
  }
  else{
    run_snv_dists <- FALSE
  }
  #checks snv_dists input
  check_snv_dists(snv_dists)
  #return the snv_dists made
  return(run_snv_dists)
}

#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_frac_intra FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
#########################################################get_clusters####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR get_clusters FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
#check that the tree object is a tree
check_tree <- function(tr){
  if(!(class(tr) == "phylo")){
    stop(paste("The tr object must be a phylogenetic tree, read into R using the ape::read.tree() function. You have supplied a ",
               class(tr)))
  }
}

#check that the pureness is a double between 1 and .5
check_pureness <- function(pureness){
  if(!(typeof(pureness) == "double")){
    stop(paste("The pureness value must be a double, you supplied type ", typeof(pureness)))
  }
  if(pureness <= 0.5 || pureness > 1){
    stop(paste("The pureness value must be between 0.5 and 1, you supplied a value of ", pureness))
  }
}


check_bootstrap <- function(bootstrap){
  #check bootstrap value is an int value or null
  if(!(is.null(bootstrap) || typeof(bootstrap) == "double")){
    stop(paste("The bootstrap value must be a double or NULL, you supplied type ", typeof(bootstrap)))
  }
  #check that the value is between 0 and 100
  if((bootstrap < 0 || bootstrap > 100) && !is.null(bootstrap)){
    stop(paste("The bootstrap value must be between 0 and 100, you supplied a value of "), bootstrap)
  }
}

#check that the names of the isolates in locs actually exist in the SNV matrix
check_tr_vs_locs <- function(tr, locs){
  #check that there are at least two in common and warn that we will subset
  if(length(intersect(names(locs), tr$tip.label)) < 2){
    stop(paste("You must supply a locs object with at least two IDs in common with your tip labels of your tr object, you have provided ",
               length(intersect(names(locs), tr$tip.label)),
               " isolates in common."))
  }
  #warn if they will be subsetting??
  if(!setequal(names(locs), tr$tip.label)){
    warning("You have not provided the same set of locs and tip labels on the tree. Will subset")
  }
}

#check get_clusters mandatory inputs
check_get_clusters_inputs <- function(tr, locs, pureness, bootstrap){
  #check that the tree is a tree
  check_tree(tr)
  #check the locs input
  check_locs(locs)
  #check tr vs. locs
  check_tr_vs_locs(tr, locs)
  #check pureness
  check_pureness(pureness)
  #check bootstrap
  check_bootstrap(bootstrap)
}

#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_clusters FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
#########################################################get_facility_fsp####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR get_facility_fsp FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
#check that the fasta file is a dna bin object
check_dna_bin <- function(fasta){
  if(all(class(fasta) != "DNAbin")){
    stop(paste("The fasta object must be of class DNAbin, you have supplied an object of class ",
               class(fasta)))
  }
}

#make sure there are at least two rownames in the fasta that match rownames in the
check_fasta_vs_locs <- function(fasta, locs){
  common_isolates <- intersect(rownames(fasta), names(locs))
  #check that there are at least 2 names of locs represented by rownames in the fasta for subsetting
  if(length(common_isolates) < 2){
    stop("Please provde a fasta object and locs object with at least two samples in common")
  }
  #check that at least two of the locs that the fasta and locs have in common appear more than once
  if(length(which(unlist(table(locs[names(locs) %in% common_isolates]) > 1))) < 2){
    stop("Please provide isolates that appear in a facility at least twice, you have not provided locs of at least two isolates in two facilities")
  }
  #if they don't have them all in common, warn that we are subsetting
  if((length(common_isolates) != length(locs)) || (length(common_isolates) != nrow(fasta))){
    warning(paste("You have provided ",
                  length(common_isolates),
                  " isolate IDs in common between locs and your fasta. Will subset."))
  }
  #warn if they are subsetting based on hospitals having less than 2 isolates
  if(any(table(locs[names(locs) %in% common_isolates]) < 2)){
    warning(paste("You have provided at least one isolate that is the only one in its location. Will subset to exclude location ",
                  paste(which(unlist(table(locs) < 2)), sep = " ", collapse = " ")))
  }
}



#do I need to check something here about how they are related? Its just confusing because how do they know what facility??
check_facility_fsp_input <- function(fasta, locs, form){
  #check form
  if(!(form == "matrix" || form == "long")){
    stop(paste("form must be either 'long' or 'matrix' to determine output format, you have provided",
               form))
  }
  #check fasta
  check_dna_bin(fasta)
  #check locs
  check_locs(locs)
  #check fasta vs locs
  check_fasta_vs_locs(fasta, locs)
}

#long form function
check_long_form_input <- function(facil_dist){
  #make sure it is a data.frame
  if(!(any(class(facil_dist) == "matrix") || any(class(facil_dist) == "data.frame"))){
    stop(paste("The fsp matrix object must be a data.frame or matrix but you provided:",
               class(facil_dist)))
  }
  #size n x n
  if(ncol(facil_dist) != nrow(facil_dist)){
    stop("The fsp matrix object must be symmetric, same number of columns as rows")
  }
  #rownames = colnames
  if(!all(rownames(facil_dist) == colnames(facil_dist))){
    stop("The fsp matrix object must be symmetric, same number of columns as rows and columns and rows must have the same names")
  }
}

#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_facility_fsp FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
##########################################################allele_freq functions####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR allele_freq functions FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
check_allele_freq_input <- function(x, subset, allele_n, alleles){
  #one option for allele_freq_between where subset !null
  #check x
  if(!(class(x) == "DNAbin" || class(x) == "raw")){
    stop(paste("The x you have provided is not a DNAbin, you provided a "),
         class(x))
  }
  #check allele_n is numeric
  if(class(allele_n) != "numeric"){
    stop(paste("The allele_n value you have provided must be numeric 1 or 2, you have provided type",
               class(allele_n)))
  }
  #check allele_n is either 1 or 2
  if(!(allele_n == 1 || allele_n == 2)){
    stop(paste("The allele_n value you have provided must be 1 or 2, you have provided",
               allele_n))
  }

  #check alleles is character vector of length two
  # & length(intersect(alleles, c("a", "c", "t", "g"))) == 2
  if(!(class(alleles) == "character" & length(alleles) == 2)){
    stop(paste("The alleles vector must be a character vector of length 2, you have provided, ",
               paste(alleles, sep = " ", collapse = " ")))
  }

  #one option for allele_freq_within where subset !null
  #check subset
  if(!is.null(subset)){
    #check subset
    if(class(subset) != "logical"){
      stop(paste("The subset vector you have provided is not logical, you provided a",
                 class(subset)))
    }
  }

}
#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR allele_freq functions FUNCTION*******************************************************#
#*******************************************************************************************************************************************#

##########################################################within_pop_var####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR within_pop_var functions FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
check_within_pop_var_inputs <- function(subset_snp_mat, subset){
  #check subset
  if(class(subset) != "logical"){
    stop(paste("The subset vector you have provided is not logical, you provided a",
               class(subset)))
  }
  #check subset_snp_mat is DNAbin
  # if(!((any(class(subset_snp_mat)) == "DNAbin" || any(class(subset_snp_mat)) == "raw") || any(class(subset_snp_mat)) == "matrixarray")){
  #   stop(paste("The subset_snp_mat you have provided is not a DNAbin, you provided a "),
  #        class(subset_snp_mat))
  # }
  #check subset_snp_mat is a matrix, DNAbin, or raw
  if(!((any(class(subset_snp_mat) == "DNAbin") ||
     any(class(subset_snp_mat) == "raw")) ||
     any(class(subset_snp_mat) == "matrixarray"))){
    stop(paste("The subset_snp_mat you have provided is not a DNAbin, you provided a "),
         class(subset_snp_mat))
  }

}
#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR within_pop_var FUNCTION*******************************************************#
#*******************************************************************************************************************************************#

#########################################################reverse_list_str####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR reverse_list_str FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
check_reverse_list_str_input <- function(ls){
  #check that it is a list
  if(class(ls) != "list"){
    stop(paste("The ls object must be a list but you provided:",
               paste(class(ls), sep = " ", collapse = " ")))
  }
  #check that the elemets of the lists are lists
  if(!all(unname(sapply(ls, class)) == "list")){
    stop(paste("The ls object must be a list of lists but you provided at lease one element that is not a list"))
  }
}
#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR reverse_list_str FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
#########################################################get_largest_subtree####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR get_largest_subtree FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
#check that subtrs are Subtrees created using ape::subtrees to look for clustering on.
#Should include all isolates of interest.
check_subtrs <- function(subtrs, isolate_labels){
  #check that it is a list
  if(class(subtrs) != "list"){
    stop(paste("The subtrs object must be the result of calling ape::subtrees, the object you have provided is not a list"))
  }
  #check that all of the elements of the list have class
  if(!all(unname(sapply(subtrs, class)) == "phylo")){
    stop(paste("The subtrs object must be the result of calling ape::subtrees, the contents of the list you provided are not of type phylo"))
  }
}

#check tip labels vs. isolate list
check_subtrs_vs_isolate_labs <- function(subtrs, isolate_labels, type){
  #check that all isolates of interest are included
  isolates <- ""
  for(i in 1:length(subtrs)){
    isolates <- c(isolates, subtrs[[i]]$tip.label)
  }
  if(length(intersect(isolates, names(isolate_labels))) != length(isolate_labels)){
    stop(paste("The subtrs object must include all of the isolates provided in the",
               type,
               "object"))
  }
}

#Named vector of labels known to cluster. Names must be equivalent to tree tip label names.
#This controls for clustering by requiring that the pure clusters must contain multiple of the control labels.
check_control_labels <- function(control_labels){
  #can be null or a named vector
  if(!is.null(control_labels) & (!is.vector(control_labels) | is.null(names(control_labels)))){
    stop(paste("control_labels must be either null (default) or a named vector of labels known to cluster"))
  }
}

check_get_largest_subtree_input <- function(subtrs, isolate_labels, control_labels, bootstrap, pureness){
  #use check_locs to check that it is a named vector of more than one isolate
  check_locs(isolate_labels)
  #check list of trees and that names are the same
  check_subtrs(subtrs)
  #check subtress vs. isolate labels
  check_subtrs_vs_isolate_labs(subtrs, isolate_labels, "isolate_labels")
  #check control labels
  check_control_labels(control_labels)
  #check control labels vs. isolate labs
  check_subtrs_vs_isolate_labs(subtrs, control_labels, "control_labels")
  #check bootstrap value
  check_bootstrap(bootstrap)
  #check pureness
  check_pureness(pureness)
}
#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_largest_subtree FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
#########################################################patient_transfer####################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR patient_transfer FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
check_thresh <- function(thresh){
  #make sure it is a single number
  if(length(thresh) != 1){
    stop(paste("thresh must be a single value, you provided",
               length(thresh)))
  }
  if(class(thresh) != "numeric"){
    stop(paste("thresh must be a numeric value, you provided",
               class(thresh)))
  }
  #make sure it is positive
  if(thresh <= 0){
    stop(paste("thresh must be a positive numeric value, you provided",
               thresh))
  }

}

#check the paths input
check_paths <- function(paths){
  #make sure it is a boolean (T or F)
  if(class(paths) != "logical"){
    stop(paste("paths argument must be a logical/logical value representing whether you want to return the shortest paths used to generate the indirect flow metric. You have provided",
               class(paths)))
  }
}


check_pt_transfer_input <- function(pt_trans_net, snv_dists, dists, locs, pt, thresh, paths){
  #check whether to run snv_dists same as get frac_intra
  #if SNV_dists doesnt exist and they didn't input locs and dists
  if(is.null(snv_dists) & is.null(dists) & is.null(locs)){
    stop("Please provide either an SNV dists matrix or dists and locs objecst so we can generate one for you.")
  }
  #make the SNV dists object if it needs to be made
  if(is.null(snv_dists)){
    #checks for get_snv_dists are done within here
    snv_dists <- get_snv_dists(dists = dists, locs = locs, pt = pt)
    run_snv_dists <- TRUE
  }
  else{
    run_snv_dists <- FALSE
  }
  #checks snv_dists input
  check_snv_dists(snv_dists)
  #check threshs input
  check_thresh(thresh)
  #if there is locs
  #check the pt_trans_net input
  check_pt_trans_net(pt_trans_net, c(snv_dists$Loc1, snv_dists$Loc2))
  #check the paths input
  check_paths(paths)

  #return the snv_dists made
  return(run_snv_dists)
}
#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR patient_transfer FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
