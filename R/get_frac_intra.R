#' Get fraction of intra-facility pairs for different snv thresholds
#'
#' @param snv_dists the output object of the get_snv_dists function
#' @param threshs SNV thresholds to use
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patient that isolate was taken from with the name being sample ID (optional)
#' @param pt_trans_net a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers' (code doesn't support missing paths, any missing paths will be represented by 0s)
#'
#' @return fraction of intra-facility pairs for different snv thresholds, lowest threshold represents lowest snv_dist in your data
#' @export
#'
#' @examples either input a snv_dists object that is the output of the get_snv_dists function or input a SNV distance matrix (made by ape::dists.dna) and a named vector of isolate locations and optionally isolate patient IDs.
get_frac_intra <- function(snv_dists = NULL, dists = NULL, locs = NULL, pt = NULL, pt_trans_net = NULL, threshs = seq(1,50,1)){

  #make one check
  run_snv_dists <- check_get_frac_intra_input(snv_dists = snv_dists, threshs = threshs, dists = dists, locs = locs, pt = pt, pt_trans_net = pt_trans_net)

  if(run_snv_dists){
    cat("Running get_snv_dists...")
    snv_dists <- get_snv_dists(dists = dists, locs = locs, pt = pt, pt_trans_net = pt_trans_net)
  }

  intra_cts <- t(sapply(threshs, function(i){
    intra <- snv_dists$Pair_Type[snv_dists$Pairwise_Dists < i]
    im <- c(i,
            table(factor(intra,levels=c('Intra-facility pair','Inter-facility pair'))),
            round(mean(intra == 'Intra-facility pair'),2),
            round(mean(intra == 'Inter-facility pair'),2))
    im
  }))
  colnames(intra_cts) <- c('Thresh','n_Intra','n_Inter','Frac_Intra','Frac_Inter')
  #remove all of the rows that have NaN as the fraction
  #intra_cts <- intra_cts %>% as.data.frame() %>% filter(!is.na(Frac_Intra), !is.na(Frac_Inter))
  intra_cts <- intra_cts %>% as.data.frame() %>% dplyr::mutate(Frac_Intra=ifelse(is.na(Frac_Intra), 0, Frac_Intra)) %>%
    dplyr::mutate(Frac_Inter=ifelse(is.na(Frac_Inter), 0, Frac_Inter))

  return(data.frame(intra_cts))
}
