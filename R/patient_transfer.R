#' Get summary of patient transfer network and closely-related isolates
#'
#' @param pt_trans_net a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers'
#' @param snv_dists the output object of the get_snv_dists function
#' @param threshs SNV thresholds to use
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patient that isolate was taken from with the name being sample ID (optional)
#'
#' @return a summary of number of closely related isolate pairs and number of patient transfers between each facility pair
#' @export
#'
#' @examples a patient transfer network and either input a snv_dists object that is the output of the get_snv_dists function or input a SNV distance matrix (made by ape::dists.dna) and a named vector of isolate locations and optionally isolate patient IDs.
patient_transfer <- function(pt_trans_net, snv_dists = NULL, dists = NULL, locs = NULL, pt = NULL, thresh = 10){
  #run checks
  run_snv_dists <- check_pt_transfer_input(pt_trans_net = pt_trans_net, snv_dists = snv_dists,
                                           dists = dists, locs = locs, pt = pt, thresh = thresh)
  #make pt_trans_net not factors
  pt_trans_net$source_facil <- as.character(pt_trans_net$source_facil)
  pt_trans_net$dest_facil <- as.character(pt_trans_net$dest_facil)
  #run get_snv_dists if necessary
  if(run_snv_dists){
    cat("Running get_snv_dists...")
    snv_dists <- get_snv_dists(dists = dists, locs = locs, pt = pt)
  }
  #subset to only facilities in the locs object
  common_locs <- intersect(unique(c(snv_dists$Loc1, snv_dists$Loc2)), unique(c(pt_trans_net$source_facil, pt_trans_net$dest_facil)))
  pt_trans_net <- pt_trans_net %>% filter(source_facil %in% common_locs,
                          dest_facil %in% common_locs)
  snv_dists <- snv_dists %>% filter(Loc1 %in% common_locs,
                                        Loc2 %in% common_locs)
  #remove any same-facility pairs
  pat_flow <- dplyr::bind_cols(pt_trans_net %>% filter(source_facil != dest_facil))
  #remove any duplicate info
  pat_flow_1 <- pat_flow[!duplicated(t(apply(pat_flow, 1, sort))),]
  #make a list for apply
  pat_flow <- split(pat_flow_1[,1:2], seq(nrow(pat_flow_1)))

  #subset to one of each pair for snv_dists
  snv_dists_sub <- subset_pairs(snv_dists)


  pat_flow_summary <- data.frame(lapply(pat_flow, function(x){
    #for each pair subset the snv_dists output to include only patients that came from those facilities
    pair <- c(as.character(x[[1]]), as.character(x[[2]]))
    #subset to only ones that are under the distance threshold
    n <- snv_dists_sub %>% filter(as.character(Loc1) %in% pair,
                                  as.character(Loc2) %in% pair,
                                  Loc1 != Loc2,
                                  as.numeric(Pairwise_Dists) <= thresh) %>%
      nrow() %>% as.numeric()
    return(n)
  }))

  #add this info back to the list
  pt_trans_summary <- data.frame(matrix(unlist(pat_flow), nrow=length(pat_flow), byrow=TRUE)) %>% mutate(facil_1 = X1, facil_2 = X2, n_closely_related_pairs = t(pat_flow_summary)) %>% select(-1, -2) %>% left_join(pat_flow_1, by = c("facil_1" = "source_facil", "facil_2" = "dest_facil"))
  pt_trans_summary$n_closely_related_pairs <- as.numeric(pt_trans_summary$n_closely_related_pairs)
  return(pt_trans_summary)
}
