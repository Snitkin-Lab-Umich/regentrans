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
  #this is where I'll do if else based on direct or indirect flow...
  if(trans_type = "indirect"){
    #run indirect_flow function, overwrite pt_trans_net
    #will have to adjust join, too
  }

  #subset to only facilities in the locs object
  common_locs <- intersect(unique(c(snv_dists$Loc1, snv_dists$Loc2)), unique(c(pt_trans_net$source_facil, pt_trans_net$dest_facil)))
  pt_trans_net <- pt_trans_net %>% dplyr::filter(source_facil %in% common_locs,
                          dest_facil %in% common_locs)
  snv_dists <- snv_dists %>% dplyr::filter(Loc1 %in% common_locs,
                                        Loc2 %in% common_locs)
  #remove any same-facility pairs
  pat_flow <- dplyr::bind_cols(pt_trans_net %>% filter(source_facil != dest_facil))

  #instead of all of that we will just left join
  pt_trans_summary <- snv_dists %>%
    mutate(Pairwise_Dists = as.numeric(Pairwise_Dists)) %>%
    #group by locations (directed)
    dplyr::group_by(Loc1, Loc2) %>%
    #summarize how many closely related isolates are in each location pair
    dplyr::summarize(n_closely_related_pairs = sum(Pairwise_Dists <= thresh)) %>%
    #remove locations that are the same
    dplyr::filter(Loc1 != Loc2) %>%
    #add in patient transfers from 1 to 2
    dplyr::left_join(pat_flow, by = c("Loc1" = "source_facil", "Loc2" = "dest_facil")) %>%
    dplyr::rename("n_1_to_2_transfers" = "n_transfers") %>%
    #add in patient transfers from 2 to 1
    dplyr::left_join(pat_flow, by = c("Loc2" = "source_facil", "Loc1" = "dest_facil")) %>%
    dplyr::rename("n_2_to_1_transfers" = "n_transfers")

  return(as.data.frame(pt_trans_summary))
}
