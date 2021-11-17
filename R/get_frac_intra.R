#' Get fraction of intra-facility pairs for different SNV thresholds
#'
#' @param pair_types the output object of the get_pair_types function
#'
#' @return fraction of intra-facility pairs for each SNV distance in dataset
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' pair_types <- get_pair_types(dists, locs)
#' get_frac_intra(pair_types)
#' }
get_frac_intra <- function(pair_types){

  #make one check
  check_get_frac_intra_input(pair_types = pair_types)

  # get intra-facility count and fraction for each pairwise SNV distance in the dataset
  intra_cts <- pair_types %>%
    dplyr::group_by(pair_type, pairwise_dist) %>%
    dplyr::mutate(pair_type = factor(pair_type, levels = c('Intra-facility pair','Inter-facility pair'))) %>%
    dplyr::count(.drop = FALSE) %>%
    tidyr::pivot_wider(names_from = pair_type, values_from = n) %>%
    dplyr::mutate(`Intra-facility pair` = ifelse('Intra-facility pair' %in% colnames(.),
                                                 as.numeric(`Intra-facility pair`), 0.0),
                  `Inter-facility pair` = ifelse('Inter-facility pair' %in% colnames(.),
                                                 as.numeric(`Inter-facility pair`), 0.0),
                  frac_intra = ifelse(`Intra-facility pair` != 0,
                                        `Intra-facility pair`/(`Intra-facility pair`+`Inter-facility pair`),
                                      0),
                  frac_inter = 1 - frac_intra
                  ) %>%
    dplyr::select(pairwise_dist, `Intra-facility pair`, `Inter-facility pair`, frac_intra, frac_inter) %>%
    dplyr::rename(n_intra = `Intra-facility pair`, n_inter = `Inter-facility pair`)
  # intra_cts

  #remove all of the rows that have NaN as the fraction
  #intra_cts <- intra_cts %>% as.data.frame() %>% dplyr::filter(!is.na(frac_intra), !is.na(frac_inter))
  intra_cts <- intra_cts %>%
    # as.data.frame() %>%
    dplyr::mutate(frac_intra=ifelse(is.na(frac_intra), 0, frac_intra)) %>%
    dplyr::mutate(frac_inter=ifelse(is.na(frac_inter), 0, frac_inter))
  return(intra_cts)
}
