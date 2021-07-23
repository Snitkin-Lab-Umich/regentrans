#' Get data frame of inter- and intra-facility pariwise SNV distances
#'
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#'
#' @return a data.frame of isolate pairs, their SNV distance, and labeled as either inter- or intra-facility pairs.
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' snv_dists <- get_snv_dists(dists, locs)
#' }
get_snv_dists <- function(dists, locs){
  #checks
  check_get_snv_dists_input(dists, locs)

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
  colnames(snps) <- c("sample1", "sample2", "pairwise_dist")

  #add locs
  snps$loc1 <- loc_sub[snps$sample1]
  snps$loc2 <- loc_sub[snps$sample2]

  snp_facility_pairs <- dplyr::bind_cols(
    snps %>% dplyr::filter(sample1 != sample2) %>%
      dplyr::mutate(pair_type=ifelse(loc1==loc2,'intra-facility pair','inter-facility pair')))

  ## should probably make this into a separate function (alphabetize loc1 and loc2)
  facil_pairs <- sapply(1:nrow(snp_facility_pairs), function(x)
    paste0(sort(c(as.character(snp_facility_pairs$loc1[x]), as.character(snp_facility_pairs$loc2[x]))), collapse = ''))

  snp_facility_pairs$loc1 <- sapply(facil_pairs, function(x) substring(x, 1, 1))
  snp_facility_pairs$loc2 <- sapply(facil_pairs, function(x) substring(x, 2, 2))
  ##

  # subset to include only one of each pair
  snp_facility_pairs <- subset_pairs(snp_facility_pairs)

  #return snp matrix
  return(snp_facility_pairs)
}


#' Subset to unique isolate pairs from a directed list
#'
#' @param snv_dists the output object of the get_snv_dists function
#'
#' @return a data.frame of isolate pairs subsetted to one row representing each pair
#' @noRd
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' snv_dists <- get_snv_dists(dists, locs)
#' subset_pairs(bind_rows(snv_dists,snv_dists))
#' }
subset_pairs <- function(snv_dists){
  #check that it is a data.frame object
  check_subset_pairs_input(snv_dists)
  #subset to one of each pair
  unique_rows <- snv_dists[!duplicated(t(apply(snv_dists, 1, sort))),]

  #return new df
  return(unique_rows)
}
