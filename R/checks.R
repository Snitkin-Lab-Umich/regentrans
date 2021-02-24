#checks for regentrans inputs

#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR get_snv_dists FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
#checks dists input to get_snv_dists function
check_dists <- function(dists){
  #if it is not a matrix
  if(!any(class(dists) == "matrix")){
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
  if(!(is.vector(locs)) | is.null(names(locs))){
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
    stop(paste("You have only supplied patient IDs for "), length(locs),
         " isolates. Please supply a named vector of patient IDs for at least 2 isolates")
  }
  #check that there are less than or equal to the number of samples in the vector than in the dists matrix
  if(length(pt) > nrow(dists)){
    stop(paste("You have supplied a list of more isolates (n = ", length(locs),
               ") with locations than exist in your SNV distance matrix (n = ",
               nrow(dists),
               ". Please make sure you have at least as many isolates in your SNV matrix as you have in your isolate location list."))
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
    warning("You have not supplied patient IDs (pt) and locations (locs) for the same samples. We will these vectors so that they contain the same isolates.")
  }
}

#check that the names of the isolates in locs actually exist in the SNV matrix
check_dists_vs_locs <- function(dists, locs){
  #check that there are less than or equal to the number of samples in the vector than in the dists matrix
  if(length(locs) > nrow(dists)){
    stop(paste("You have supplied a list of more isolates (n = ", length(locs),
               ") with locations than exist in your SNV distance matrix (n = ",
               nrow(dists),
               ". Please make sure you have at least as many isolates in your SNV matrix as you have in your isolate location list."))
  }
  #check if the names of the locs isolates are a subset of the names of the dist matrix isolates
  if(!all(names(locs) %in% rownames(dists))){
    stop("Some of the isolates you have provided locations for are not in the SNV distance matrix (dists). Please subset the locs vector to include only isolates in the SNV distance matrix (dists).")
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

#wrapper for if pt is included
check_get_snv_dists_input_pt <- function(dists, locs, pt){
  #check that the pt object is a named vector
  check_pt(pt, dists)
  #check that the pt and dists objects have the same lengths
  check_pt_vs_locs(pt, locs)
}

#wrapper function if pt not included
check_get_snv_dists_input_no_pt <- function(dists, locs){
  #check that the dists object is the snv object returned by dist.dna
  check_dists(dists)
  #check that the locs object is a named vector
  check_locs(locs)
  #check that the locs names exist in the dists dataframe
  check_dists_vs_locs(dists, locs)
}

#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_snv_dists FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
#############################################################################################################################################
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
  if(!(ncol(snv_dists) == 8 | ncol(snv_dists) == 6)){
    stop(paste("The snv_dists object must be the output of the get_snv_dists() function, but you provided a data.frame with ",
               ncol(snv_dists), " columns."))
  }
  #check the colnames
  if(!((ncol(snv_dists) == 8 &&
        all(colnames(snv_dists) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2", "Patient1",  "Patient2", "Pair_Type"))) ||
       (ncol(snv_dists) == 6 &&
        all(colnames(snv_dists) == c("Isolate1", "Isolate2", "Pairwise_Dists", "Loc1", "Loc2" "Pair_Type"))))){
    stop(paste("The snv_dists object must be the output of the get_snv_dists() function, but the data.frame you provided has ",
               ncol(snv_dists), " columns that are not the output columns needed."))
  }
  #check that there are at least one row
  if(nrow(snv_dists) < 1){
      stop(paste("Your snv_dists input has ", nrow(snv_dists), " facility pairs. Please use an snv_dists input that has 1 or more pairs (rows)"))
  }
  #check pairwise dist column is numeric
  if(!class(snv_dists$Pairwise_Dists) == "numeric"){
    stop(paste("Your snv_dists input does not have numeric pairwise distances, you supplied one with type",
               class(snv_dists$Pairwise_Dists)))
  }

}

#some vector of numbers, max number isn't > max snv distance or negative
check_threshs <- function(threshs, snv_dists){
  #check numeric vector
  if(!(is.vector(threshs) && class(threshs) == "numeric")){
    stop(paste("threshs must be a numeric vector, you provided a ", typeof(threshs)))
  }
  #check that max number isn't > max snv distance or negative
  #we are not meeting this condition right now so I will just throw a warning
  if(!max(threshs) <= max(snv_dists$Pairwise_dists)){
    stop(paste("Your max threshold ", max(threshs), " is greater than your max CNV distance of ", max(snv_dists$Pairwise_dists),
                  "If you are using the default thereshold input, please enter a numeric vector of thresholds"))
  }
  if(any(threshs < 0)){
    stop("You provided a threshold below 0")
  }
}

check_get_frac_intra_input <- function(snv_dists, threshs){
  #checks snv_dists input
  check_snv_dists(snv_dists)
  #check threshs input
  check_threshs(threshs, snv_dists)
}

#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_frac_intra FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
#############################################################################################################################################
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
  if(bootstrap < 0 || bootstrap > 100){
    stop(paste("The bootstrap value must be between 0 and 100, you supplied a value of"), bootstrap)
  }
}

#check that the names of the isolates in locs actually exist in the SNV matrix
check_tr_vs_locs <- function(tr, locs){
  #check that there are less than or equal to the number of samples in the vector than in the tree tip labels
  if(length(locs) > length(tr$tip.label)){
    stop(paste("You have supplied a list of more isolates (n = ", length(locs),
               ") with locations than exist in your tree (n = ",
               length(tr$tip.label),
               ". Please make sure you have at least as many isolates in your SNV matrix as you have in your isolate location list."))
  }
  #check if the names of the locs isolates are a subset of the names of the dist matrix isolates
  if(!all(names(locs) %in% tr$tip.label)){
    stop("Some of the isolates you have provided locations for are not in the SNV distance matrix (dists). Please subset the locs vector to include only isolates in the SNV distance matrix (dists).")
  }
  #warn if they will be subsetting??
  if(!setequal(names(locs), tr$tip.label)){
    warning("You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset")
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
#############################################################################################################################################
#############################################################################################################################################
#*******************************************************************************************************************************************#
#***********************************************START CHECKS FOR get_facility_fsp FUNCTION*****************************************************#
#*******************************************************************************************************************************************#
#check that the fasta file is a dna bin object
check_dna_bin <- function(fasta){
  if(class(fasta) != "DNAbin"){
    stop(paste("The fasta object must be of class DNAbin, you have supplied an object of class ",
               class(fasta)))
  }
}

#make sure there are at least two rownames in the fasta that match rownames in the
check_fasta_vs_locs <- function(fasta, locs){
  #check that there are at least 2 names of locs represented by rownames in the fasta for subsetting
  if(length(intersect(rownames(fasta), names(locs))) < 2){
    stop("Please provde a fasta object and locs object with at least two samples in common")
  }
}

#do I need to check something here about how they are related? Its just confusing because how do they know what facility??
check_facility_fsp <- function(fasta, locs){
  #check fasta
  check_dna_bin(fasta)
  #check locs
  check_locs(locs)
  #check fasta vs locs
  check_fasta_vs_locs(fasta, locs)
}


#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_facility_fsp FUNCTION*******************************************************#
#*******************************************************************************************************************************************#
#############################################################################################################################################

