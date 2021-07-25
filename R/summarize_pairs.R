#' Summarize extent of relatedness of intra- and inter-facility pairs by facility pair
#'
#' @param snv_dists the output object of the get_snv_dists function
#' @param summary_fns vector of summary functions for pairwise distances as character strings (default: c("min"))
#' @param threshs SNV thresholds to use for pairwise distances (default: seq(5, 20, 5))
#'
#' @return a summary of isolate pairs between each facility pair
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' snv_dists <- get_snv_dists(dists, locs)
#' summarize_pairs(snv_dists = snv_dists)
#' }
summarize_pairs <- function(snv_dists, summary_fns = c("min"), threshs = seq(5,20,5)){
  #run checks
  check_summarize_pairs_input(snv_dists = snv_dists, summary_fns = summary_fns, threshs = threshs)

  inter_dists_stats_summary <- data.frame(loc1=character(), loc2=character())

  if(!is.null(summary_fns)){
    inter_dists_stats_summary <- lapply(summary_fns, function(x){
      snv_dists %>% dplyr::group_by(loc1, loc2) %>%
        dplyr::summarize(!!dplyr::quo_name(paste0('dist_', x)) := get(x)(pairwise_dist), .groups = 'keep')
    }) %>%
      purrr::reduce(dplyr::full_join, by = c("loc1", "loc2"))
  }

  inter_dists_thresh_summary <- data.frame(loc1=character(), loc2=character())
  if(!is.null(threshs)){
    inter_dists_thresh_summary <- lapply(threshs, function(x){
      snv_dists %>% dplyr::group_by(loc1, loc2) %>%
        dplyr::summarize(!!dplyr::quo_name(paste0('under_', x)) := sum(pairwise_dist < x), .groups = 'keep')
    }) %>%
      purrr::reduce(dplyr::full_join, by = c("loc1", "loc2"))
  }

  inter_dists_summary <- dplyr::full_join(inter_dists_stats_summary, inter_dists_thresh_summary,
                                     by = c("loc1", "loc2"))

  return(inter_dists_summary)
}


#' Merge summarized data about facility pairs
#'
#' @param patient_flow output of get_patient_flow function
#' @param isolate_pair_summary output of summarize_pairs function
#' @param fsp_long long-form output of get_facility_fsp function
#'
#' @return merged dataframe
#' @export
#' @description Merge summarized data about facility pairs including patient flow, isolate pairs, and gene flow (Fsp).
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' snv_dists <- get_snv_dists(dists, locs)
#' isolate_pair_summary <- summarize_pairs(snv_dists)
#' patient_flow <- get_patient_flow(edge_df = pt_trans_df)
#' fsp_long <- make_long_form(fsp)
#' merge_inter_summaries(patient_flow, isolate_pair_summary, fsp_long)
#' }
merge_inter_summaries <- function(patient_flow = NULL, isolate_pair_summary = NULL, fsp_long = NULL){
  check_merge_inter_summaries_input(patient_flow = patient_flow, isolate_pair_summary = isolate_pair_summary, fsp_long = fsp_long)

  if(is.null(patient_flow)) patient_flow <- data.frame(loc1=character(), loc2=character())
  if(is.null(isolate_pair_summary)) isolate_pair_summary <- data.frame(loc1=character(), loc2=character())
  if(is.null(fsp_long)) fsp_long <- data.frame(loc1=character(), loc2=character())

  merged_df <- dplyr::full_join(isolate_pair_summary, fsp_long, by = c("loc1", "loc2")) %>%
    dplyr::left_join(patient_flow, by = c("loc1", "loc2"))

  return(merged_df)

}
