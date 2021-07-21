#' Get number of closely related pairs between facilities for different SNV distance thresholds
#'
#' @param snv_dists the output object of the get_snv_dists function
#' @param summary_fns vector of summary functions for pairwise distances as character strings (default: c("min"))
#' @param threshs SNV thresholds to use for pairwise distances (default: seq(5, 20, 5))
#'
#' @return a summary of isolate pairs between each facility pair.
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' snv_dists <- get_snv_dists(dists, locs)
#' summarize_inter_pairs(snv_dists = snv_dists)
#' }
summarize_inter_pairs <- function(snv_dists, summary_fns = c("min"), threshs = seq(5,20,5)){
  #run checks
  check_summarize_inter_pairs_input(snv_dists = snv_dists, summary_fns = summary_fns, threshs = threshs)

  val_funs <- lapply(threshs, function(x) (function(a) as.integer(a < x)))
  names(val_funs) <- paste0("<", threshs)

  inter_dists_stats_summary <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Loc1", "Loc2"))

  if(!is.null(summary_fns)){
    inter_dists_stats_summary <- lapply(summary_fns, function(x){
      snv_dists %>% group_by(Loc1, Loc2) %>%
        summarize(!!quo_name(paste0('dists_', x)) := get(x)(Pairwise_Dists), .groups = 'keep')
    }) %>%
      purrr::reduce(full_join, by = c("Loc1", "Loc2"))
  }

  inter_dists_thresh_summary <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Loc1", "Loc2"))
  if(!is.null(threshs)){
    inter_dists_thresh_summary <- lapply(threshs, function(x){
      snv_dists %>% group_by(Loc1, Loc2) %>%
        summarize(!!quo_name(paste0('under_', x)) := sum(Pairwise_Dists < x), .groups = 'keep')
    }) %>%
      purrr::reduce(full_join, by = c("Loc1", "Loc2"))
  }

  inter_dists_summary <- full_join(inter_dists_stats_summary, inter_dists_thresh_summary,
                                     by = c("Loc1", "Loc2"))

  return(inter_dists_summary)
}


#' Merge summarized data about facility pairs
#'
#' @param patient_flow output of get_patient_flow function
#' @param inter_pair_summary output of summarize_inter_pairs function
#' @param fsp output of get_facility_fsp function
#'
#' @return merged dataframe
#' @export
#'
#' @examples
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' snv_dists <- get_snv_dists(dists, locs)
#' inter_pair_summary <- summarize_inter_pairs(snv_dists)
#' patient_flow <- get_patient_flow(edge_df = pt_trans_df)
#' merge_inter_summaries(patient_flow, inter_pair_summary, fsp)
merge_inter_summaries <- function(patient_flow = NULL, inter_pair_summary = NULL, fsp = NULL){
  check_merge_inter_summaries_input(patient_flow = patient_flow, inter_pair_summary = inter_pair_summary, fsp = fsp)

  if(is.null(patient_flow)) patient_flow <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Loc1", "Loc2"))
  if(is.null(inter_pair_summary)) inter_pair_summary <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Loc1", "Loc2"))
  if(is.null(fsp)) fsp <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Loc1", "Loc2"))

  merged_df <- full_join(inter_pair_summary, fsp) %>% left_join(patient_flow)

  return(merged_df)

}
