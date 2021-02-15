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

}

#checks locs input to get_snv_dists function
check_locs <- function(locs, dists){
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
  #check that there are less than or equal to the number of samples in the vector than in the dists matrix
  if(length(locs) > nrow(dists)){
    stop(paste("You have supplied a list of more isolates (n = ", length(locs), 
               ") with locations than exist in your SNV distance matrix (n = ", 
               nrow(dists), 
               ". Please make sure you have at least as many isolates in your SNV matrix as you have in your isolate location list."))
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
  #here would I want to subet it myself
  if(length(pt) != length(locs)){
    stop(paste("You have supplied a patient vector of length ", length(pt), " and  location vector of length ", length(locs), ". Please subset these lists so that they have the same number of isolates."))
  }
  #check that the names match between pt and locs
  if(!setequal(names(pt), names(locs))){
    stop("You have not supplied patient IDs (pt) and locations (locs) for the same samples. Pleasesubset these vectors so that they contain the same isolates.")
  }
}

#check that the names of the isolates in locs actually exist in the SNV matrix
check_dists_vs_locs <- function(dists, locs){
  #check if the names of the locs isolates are a subset of the names of the dist matrix isolates
  if(!all(names(locs) %in% rownames(dists))){
    stop("Some of the isolates you have provided locations for are not in the SNV distance matrix (dists). Please subset the locs vector to include only isolates in the SNV distance matrix (dists).")
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
  check_locs(locs, dists)
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



#*******************************************************************************************************************************************#
#***********************************************END CHECKS FOR get_frac_intra FUNCTION*******************************************************#
#*******************************************************************************************************************************************#

