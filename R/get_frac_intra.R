#' Get fraction of intra-facility pairs for different SNV thresholds
#'
#' @param snv_dists the output object of the get_snv_dists function
#'
#' @return fraction of intra-facility pairs for each SNV distance in dataset
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(sample_id, facility) %>% tibble::deframe()
#' snv_dists <- get_snv_dists(dists, locs)
#' get_frac_intra(snv_dists)
#' }
get_frac_intra <- function(snv_dists){

  #make one check
  check_get_frac_intra_input(snv_dists = snv_dists)

  # get intra-facility count and fraction for each pairwise SNV distance in the dataset
  intra_cts <- snv_dists %>%
    dplyr::group_by(pair_type, pairwise_dist) %>%
    dplyr::mutate(pair_type = factor(pair_type, levels = c('intra-facility pair','inter-facility pair'))) %>%
    dplyr::count(.drop = FALSE) %>%
    tidyr::pivot_wider(names_from = pair_type, values_from = n) %>%
    dplyr::mutate(`intra-facility pair` = ifelse('intra-facility pair' %in% colnames(.),
                                                 as.numeric(`intra-facility pair`), 0.0),
                  `inter-facility pair` = ifelse('inter-facility pair' %in% colnames(.),
                                                 as.numeric(`inter-facility pair`), 0.0),
                  frac_intra = ifelse(`intra-facility pair` != 0,
                                        `intra-facility pair`/(`intra-facility pair`+`inter-facility pair`),
                                      0),
                  frac_inter = 1 - frac_intra
                  ) %>%
    dplyr::select(pairwise_dist, `intra-facility pair`, `inter-facility pair`, frac_intra, frac_inter) %>%
    dplyr::rename(n_intra = `intra-facility pair`, n_inter = `inter-facility pair`)
  # intra_cts

  #remove all of the rows that have NaN as the fraction
  #intra_cts <- intra_cts %>% as.data.frame() %>% dplyr::filter(!is.na(frac_intra), !is.na(frac_inter))
  intra_cts <- intra_cts %>%
    # as.data.frame() %>%
    dplyr::mutate(frac_intra=ifelse(is.na(frac_intra), 0, frac_intra)) %>%
    dplyr::mutate(frac_inter=ifelse(is.na(frac_inter), 0, frac_inter))
  return(intra_cts)
}
