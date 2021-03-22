
#' Get summary of patient transfer network and closely-related isolates
#'
#' @param pt_trans_net a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers'
#' @param snv_dists the output object of the get_snv_dists function
#' @param threshs SNV thresholds to use
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patient that isolate was taken from with the name being sample ID (optional)
#'
#' @return fraction of intra-facility pairs for different snv thresholds, lowest threshold represents lowest snv_dist in your data
#' @export
#'
#' @examples a patient transfer network and either input a snv_dists object that is the output of the get_snv_dists function or input a SNV distance matrix (made by ape::dists.dna) and a named vector of isolate locations and optionally isolate patient IDs.
patient_transfer <- function(pt_trans_net, snv_dists = NULL, dists = NULL, locs = NULL, pt = NULL, thresh = 10){
  #run checks
  check_pt_transfer_input(pt_trans_net, snv_dists, dists, locs, pt, thresh)
  #run get_snv_dists if necessary

  pat_flow <- dplyr::bind_cols(pt_trans_net %>% filter(source_facil != dest_facil))
  pat_flow_1 <- pat_flow[!duplicated(t(apply(pat_flow, 1, sort))),]
  pat_flow <- split(pat_flow_1[,1:2], seq(nrow(pat_flow_1)))

  #subset to one of each pair
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
  df <- data.frame(matrix(unlist(pat_flow), nrow=length(pat_flow), byrow=TRUE)) %>% mutate(facil_1 = X1, facil_2 = X2, n_closely_related_pairs = t(pat_flow_summary)) %>% select(-1, -2) %>% left_join(pat_flow_1, by = c("facil_1" = "source_facil", "facil_2" = "dest_facil"))

}

#make a source destination pair test matrix
# mat <- data.frame(matrix(data = c(0, 20, 30,
#                         20, 0, 26,
#                         30, 26, 0), nrow = 3, ncol = 3))
# rownames(mat) <- c("A", "B", "C")
# colnames(mat) <- c("A", "B", "C")
# pat_flow <- na.omit(data.frame(as.table(as.matrix(mat))))
# #pat_flow <- dplyr::bind_cols(pat_flow %>% filter(Var1 != Var2))
# colnames(pat_flow) <- c("source_facil", "dest_facil", "n_transfers")
# pat_flow$n_transfers <- as.numeric(pat_flow$n_transfers)
#
# #try to do stuff with this now to snv_dists output
# snv_dists %>% dplyr::left_join(pat_flow, by = c("Loc1" = "source_facil", "Loc2" = "dest_facil"))
#
#
# #write function to summarize
# #facil_1, Facil_2, transfer_info, n_closely_related_pairs
# #input: snv_dists output (or stuff to make it), patient flow network, threshold
# thresh = 50
#
# #subset to the one of each pair
# pat_flow <- dplyr::bind_cols(pat_flow %>% filter(source_facil != dest_facil))
# pat_flow_1 <- pat_flow[!duplicated(t(apply(pat_flow, 1, sort))),]
# pat_flow <- split(pat_flow_1[,1:2], seq(nrow(pat_flow_1)))
#
# #subset to one of each pair
# #snv_dists_sub <- subset_pairs(snv_dists)
#
# pat_flow_summary <- data.frame(lapply(pat_flow, function(x){
#   #for each pair subset the snv_dists output to include only patients that came from those facilities
#   pair <- c(as.character(x[[1]]), as.character(x[[2]]))
#   n <- snv_dists_sub %>% filter(as.character(Loc1) %in% pair,
#                                 as.character(Loc2) %in% pair,
#                                 Loc1 != Loc2,
#                                 as.numeric(Pairwise_Dists) <= thresh) %>%
#     nrow() %>% as.numeric()
#   return(n)
#   #subset to only ones that are under the distance threshold
#
# }))
#
# #add this info back to the list
# df <- data.frame(matrix(unlist(pat_flow), nrow=length(pat_flow), byrow=TRUE)) %>% mutate(facil_1 = X1, facil_2 = X2, n_closely_related_pairs = t(pat_flow_summary)) %>% select(-1, -2) %>% left_join(pat_flow_1, by = c("facil_1" = "source_facil", "facil_2" = "dest_facil"))
