# get fraction of intra-facility pairs for different snv thresholds
# dists - output of get_snv_dists - will want to change either the name of this argument or the other one so the arguments always mean the same input
# threshs - what snv thresholds to use #some vector of numbers, max number isn't > max snv distance or negative
# what if we call it snv_dists? because the other one is getting snv_dists?

#' Get fraction of intra-facility pairs for different snv thresholds
#'
#' @param snv_dists the output object of the get_snv_dists function
#' @param threshs SNV thresholds to use
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patient that isolate was taken from with the name being sample ID (optional)
#'
#' @return fraction of intra-facility pairs for different snv thresholds
#' @export
#'
#' @examples either input a snv_dists object that is the output of the get_snv_dists function or input a SNV distance matrix (made by ape::dists.dna) and a named vector of isolate locations and optionally isolate patient IDs.
get_frac_intra <- function(snv_dists = NULL, dists = NULL, locs = NULL, pt = NULL, threshs = seq(1,50,1)){

  #if there is a snv_dists input, check inputs
  if(!is.null(snv_dists)){
    check_get_frac_intra_input(snv_dists, threshs)
  }
  else{
    #if there isn't run get_snv_dists
    snv_dists <- get_snv_dists(dists, locs, pt)
    check_threshs(threshs, snv_dists)
  }

  intra_cts <- t(sapply(threshs, function(i){
    intra <- snv_dists$intra[snv_dists$Pairwise_dists < i & !is.na(snv_dists$intra)]
    im <- c(i,
            table(factor(intra,levels=c('Intra-facility pair','Inter-facility pair'))),
            round(mean(intra == 'Intra-facility pair'),2),
            round(mean(intra == 'Inter-facility pair'),2))
    names(im) <- c('thresh','inter','intra','frac_intra','frac_inter')
    im
  }))
  return(data.frame(intra_cts))
}

#might want to change the name/call the get_snv_dists
#within it just so we don't have to run them in order,
#could write a wrapper function that calls both - get_frac_intra_from_snv_dists
