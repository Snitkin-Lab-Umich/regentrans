#' Get fraction of intra-facility pairs for different snv thresholds
#'
#' @param snv_dists the output object of the get_snv_dists function
#'
#' @return fraction of intra-facility pairs for different snv thresholds, lowest threshold represents lowest snv_dist in your data
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' get_frac_intra(dists = dists, locs = locs)
#' }
get_frac_intra <- function(snv_dists){

  #make one check
  check_get_frac_intra_input(snv_dists = snv_dists)

  # get intra-facility count and fraction for each pairwise SNV distance in the dataset
  intra_cts <- snv_dists %>%
    dplyr::group_by(Pair_Type, Pairwise_Dists) %>%
    dplyr::mutate(Pair_Type = factor(Pair_Type, levels = c('Intra-facility pair','Inter-facility pair'))) %>%
    dplyr::count(.drop = FALSE) %>%
    tidyr::pivot_wider(names_from = Pair_Type, values_from = n) %>%
    dplyr::mutate(`Intra-facility pair` = ifelse('Intra-facility pair' %in% colnames(.),
                                                 as.numeric(`Intra-facility pair`), 0.0),
                  `Inter-facility pair` = ifelse('Inter-facility pair' %in% colnames(.),
                                                 as.numeric(`Inter-facility pair`), 0.0),
                  Frac_Intra = ifelse(`Intra-facility pair` != 0,
                                        `Intra-facility pair`/(`Intra-facility pair`+`Inter-facility pair`),
                                      0),
                  Frac_Inter = 1 - Frac_Intra
                  ) %>%
    dplyr::select(Pairwise_Dists, `Intra-facility pair`, `Inter-facility pair`, Frac_Intra, Frac_Inter) %>%
    dplyr::rename(n_Intra = `Intra-facility pair`, n_Inter = `Inter-facility pair`)
  # intra_cts

  #remove all of the rows that have NaN as the fraction
  #intra_cts <- intra_cts %>% as.data.frame() %>% dplyr::filter(!is.na(Frac_Intra), !is.na(Frac_Inter))
  intra_cts <- intra_cts %>%
    # as.data.frame() %>%
    dplyr::mutate(Frac_Intra=ifelse(is.na(Frac_Intra), 0, Frac_Intra)) %>%
    dplyr::mutate(Frac_Inter=ifelse(is.na(Frac_Inter), 0, Frac_Inter))
  return(intra_cts)
}
