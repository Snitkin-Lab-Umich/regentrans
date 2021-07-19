#' Get summary of patient transfer network and closely-related isolates
#'
#' @param pt_trans_net a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers' (code doesn't support missing paths, any missing paths will be represented by 0s)
#' @param snv_dists the output object of the get_snv_dists function
#' @param thresh SNV thresholds to use
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patient that isolate was taken from with the name being sample ID (optional)
#' @param paths boolean value, TRUE if you want the shortest paths returned, FALSE if you don't
#'
#' @return a summary of number of closely related isolate pairs and number of direct patient transfers and indirect flow metrics between each facility pair. If paths = TRUE, a list of summary (pt_trans_summary) and shortest paths used (paths).
#' @export
#'
#' @examples get_patient_transfers(pt_trans_net = pt_flow, dists = dists, locs = locs)
get_patient_transfers <- function(pt_trans_net, snv_dists = NULL, dists = NULL, locs = NULL, pt = NULL, thresh = 10, paths = FALSE){
  #run checks
  run_snv_dists <- check_pt_transfer_input(pt_trans_net = pt_trans_net, snv_dists = snv_dists,
                                           dists = dists, locs = locs, pt = pt, thresh = thresh, paths = paths)

  #make pt_trans_net not factors
  pt_trans_net$source_facil <- as.character(pt_trans_net$source_facil)
  pt_trans_net$dest_facil <- as.character(pt_trans_net$dest_facil)

  # fill in missing source and destination facilities (doesn't change results, but will error out otherwise)
  all_facils <- unique(c(pt_trans_net$source_facil,pt_trans_net$dest_facil))
  not_in_source <- all_facils[!(all_facils %in% pt_trans_net$source_facil)]
  not_in_dest <- all_facils[!(all_facils %in% pt_trans_net$dest_facil)]
  if(length(not_in_source) != 0 | length(not_in_dest) != 0){
    pt_trans_net <- dplyr::bind_rows(pt_trans_net,
                                     dplyr::bind_cols(source_facil = not_in_source,
                                                      dest_facil = pt_trans_net$dest_facil[1],
                                                      n_transfers = NA),
                                     dplyr::bind_cols(source_facil = pt_trans_net$source_facil[1],
                                                      dest_facil = not_in_dest,
                                                      n_transfers = NA))
    pt_trans_net <- pt_trans_net %>% tidyr::expand(source_facil, dest_facil) %>%
      dplyr::left_join(pt_trans_net, by = c("source_facil", "dest_facil"))
  }

  #run get_snv_dists if necessary
  if(run_snv_dists){
    message("Running get_snv_dists...")
    snv_dists <- get_snv_dists(dists = dists, locs = locs, pt = pt)
  }

  #run indirect flow
  ind_flow_output <- get_indirect_flow(pt_trans_net)
  pt_trans_net_i <- ind_flow_output$transfer_network
  paths_list <- ind_flow_output$paths

  #subset to only facilities in the locs object
  common_locs <- intersect(unique(c(snv_dists$Loc1, snv_dists$Loc2)), unique(c(pt_trans_net$source_facil, pt_trans_net$dest_facil)))
  pt_trans_net <- pt_trans_net %>% dplyr::filter(source_facil %in% common_locs,
                          dest_facil %in% common_locs)
  snv_dists <- snv_dists %>% dplyr::filter(Loc1 %in% common_locs,
                                        Loc2 %in% common_locs)

  #remove any same-facility pairs
  pat_flow <- dplyr::bind_cols(pt_trans_net %>% dplyr::filter(source_facil != dest_facil))

  #instead of all of that we will just left join
  pt_trans_summary <- snv_dists %>%
    dplyr::mutate(Pairwise_Dists = as.numeric(Pairwise_Dists)) %>%
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
    dplyr::rename("n_2_to_1_transfers" = "n_transfers") %>%
    #add in indirect flow 1 to 2
    dplyr::left_join(pt_trans_net_i, by = c("Loc1" = "source_facil", "Loc2" = "dest_facil")) %>%
    dplyr::rename("indirect_flow_metric_1_to_2" = "pt_trans_metric") %>%
    #add in indirect flow 2 to 1
    dplyr::left_join(pt_trans_net_i, by = c("Loc2" = "source_facil", "Loc1" = "dest_facil")) %>%
    dplyr::rename("indirect_flow_metric_2_to_1" = "pt_trans_metric")

  #make df
  pt_trans_summary <- as.data.frame(pt_trans_summary)

  if(paths == TRUE){
    #return paths and summary as a list
    returns <- list("pt_trans_summary" = pt_trans_summary, "paths" = paths_list)

    return(returns)
  }
  else{
    return(pt_trans_summary)
  }
}


#' Calculate indirect patient flow from patient transfer network
#'
#' @param pt_trans_net a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers' (code doesn't support missing paths, any missing paths will be represented by 0s)
#'
#' @return facility x facility matrix of metric of patient flow between each facility pair
#' @noRd
#'
#' @examples get_indirect_flow(pt_trans_net)
get_indirect_flow <- function(pt_trans_net){
  #don't want to subset before getting here, need whole network for indirect
  #checks
  check_pt_trans_net(pt_trans_net, unique(c(as.character(pt_trans_net$source_facil), as.character(pt_trans_net$dest_facil))))

  #make matrix format
  trans_mat <- tidyr::pivot_wider(pt_trans_net, names_from = source_facil, values_from = n_transfers)
  trans_mat <- as.data.frame(trans_mat[,2:ncol(trans_mat)])
  rownames(trans_mat) <- colnames(trans_mat)
  trans_mat = t(trans_mat)

  #make graph
  g <- igraph::graph_from_adjacency_matrix(as.matrix(trans_mat),mode='directed',weighted = TRUE)

  #name nodes in network that we have data for

  #modify edge weights from n facilities -> normalize/invert
  out_strength = igraph::strength(g,mode='out') # get number of outgoing patient transfers for each vertex
  tail_vert = igraph::tail_of(g,igraph::E(g)) # get tail (source) vertex for each edge
  edwt_sum = sapply(names(tail_vert), function(x) out_strength[names(out_strength) == x]) # get number of outgoing patient transfers of tail vertex for each edge
  igraph::E(g)$weight = -log10(igraph::E(g)$weight/edwt_sum) # normalize edge weight by number of outgoing patient transfers of source vertex and take negative log (to use to calculate shortest paths)

  #maybe subset to nodes that have info in our dataset

  #find shortest path function -> igraph::shortest.paths()
  sp <- g %>% igraph::shortest.paths(mode="out") %>% as.data.frame()

  #make long form
  trans_net_i <- sp %>% tibble::as_tibble() %>% dplyr::mutate(source_facil = colnames(sp)) %>% tidyr::pivot_longer(!source_facil, names_to = "dest_facil", values_to = "pt_trans_metric")

  #make them each -(10^x)
  trans_net_i$pt_trans_metric <- 10^(-trans_net_i$pt_trans_metric)
  trans_net_i$pt_trans_metric[is.na(trans_net_i$pt_trans_metric)] <- 0

  #if they asked for the paths, find them and make/return list
  paths <- igraph::get.shortest.paths(g, 1, mode = "out")
  returns <- list("transfer_network" = trans_net_i, "paths" = paths)
  return(returns)
}

