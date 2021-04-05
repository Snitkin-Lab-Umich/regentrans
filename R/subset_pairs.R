#write subset_pairs



#' subset_pairs
#'
#' @param snv_dists the output object of the get_snv_dists function
#'
#' @return a data.frame of isolate pairs subsetted to one row representing each pair
#' @export
#'
#' @examples subset_pairs(snv_dists)
subset_pairs <- function(snv_dists){
  #check that it is a data.frame object
  #check_snv_dists(snv_dists)
  #subset to one of each pair
  unique_rows <- snv_dists[!duplicated(t(apply(snv_dists, 1, sort))),]
  #return new df
  return(unique_rows)
}
