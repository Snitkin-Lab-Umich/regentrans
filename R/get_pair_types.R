#' Get data frame of inter- and intra-facility pariwise SNV distances
#'
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patients each isolate originated from, with the name being the sample ID. If this information is unavailable, set pt = NULL.
#'
#' @return a data.frame of isolate pairs, their SNV distance, and labeled as either inter- or intra-facility pairs.
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' pt <- metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()
#' pair_types <- get_pair_types(dists, locs, pt)
#' }
get_pair_types <- function(dists, locs, pt){
  #checks
  #TO DO: add a check for pt (either null or a named vector of isolates)
  check_get_pair_types_input(dists, locs)

  #make the subsetted isolates object
  isolates <- intersect(names(locs), rownames(dists))

  #subset by locs
  #list ones in common before subsetting
  loc_sub <- locs[isolates]

  #subset dists to isolates
  dists_sub <- dists[isolates, isolates]

  #make df
  snps <- stats::na.omit(data.frame(as.table(as.matrix(dists_sub))))
  #change freq colname?
  colnames(snps) <- c("isolate1", "isolate2", "pairwise_dist")

  #add locs
  snps$loc1 <- loc_sub[snps$isolate1]
  snps$loc2 <- loc_sub[snps$isolate2]

  #add pts
  if(!is.null(pt)){
    snps$pt1 <- pt[snps$isolate1]
    snps$pt2 <- pt[snps$isolate2]
  }

  snp_facility_pairs <- dplyr::bind_cols(
    snps %>% dplyr::filter(isolate1 != isolate2) %>%
      dplyr::mutate(pair_type=ifelse(loc1==loc2,'Intra-facility pair','Inter-facility pair')))

  # subset to include only one of each pair
  snp_facility_pairs <- snp_facility_pairs %>% dplyr::arrange(pt1) %>% subset_pairs()

  #if there is pt info
  if(!is.null(pt)){
    snp_facility_pairs <- snp_facility_pairs %>% group_by(loc1, loc2, pt1, pt2) %>% slice(which.min(pairwise_dist))
  }

  #return snp matrix
  return(snp_facility_pairs %>% ungroup() %>% as.data.frame())
}


#' Subset to unique isolate pairs from a directed list
#'
#' @param pair_types the output object of the get_pair_types function
#'
#' @return a data.frame of isolate pairs subsetted to one row representing each pair
#' @noRd
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' pair_types <- get_pair_types(dists, locs)
#' subset_pairs(bind_rows(pair_types,pair_types))
#' }
subset_pairs <- function(pair_types){
  #check that it is a data.frame object
  check_subset_pairs_input(pair_types)
  #subset to one of each pair
  unique_rows <- pair_types[!duplicated(t(apply(pair_types, 1, sort))),]

  #return new df
  return(unique_rows)
}
