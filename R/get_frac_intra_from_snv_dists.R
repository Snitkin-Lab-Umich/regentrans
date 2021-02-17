#wrapper function for get_frac_intra and get_snp_dists functions
#would we want to return both outputs in a lsit or is the frac_intra the thing of interest here?

#' Get fraction of intra-facility pairs for different snv thresholds from SNV distance matrix
#'
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patient that isolate was taken from with the name being sample ID (optional)
#' @param threshs SNV thresholds to use
#'
#' @return Fraction of intra-facility pairs for different snv thresholds
#' @export
#'
#' @examples
get_frac_intra_from_snv_dists <- function(dists, locs, pt, threshs = seq(1,50,1)){
  #run get_frac_intra
  snv_dists <- get_snv_dists(dists, locs, pt)
  frac_intra <- get_frac_intra(snv_dists, threshs)
  return(frac_intra)
}
