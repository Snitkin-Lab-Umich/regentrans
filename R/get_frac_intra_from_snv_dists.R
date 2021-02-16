#wrapper function for get_frac_intra and get_snp_dists functions
#would we want to return both outputs in a lsit or is the frac_intra the thing of interest here? 
get_frac_intra_from_snv_dists <- function(dists, locs, pt, threshs = seq(1,50,1)){
  #run get_frac_intra
  snv_dists <- get_snv_dists(dists, locs, pt)
  frac_intra <- get_frac_intra(snv_dists, threshs)
  return(frac_intra)
}