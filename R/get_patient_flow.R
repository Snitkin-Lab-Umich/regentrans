
#' Summarize inter-facility patient transfer network
#'
#' @param edge_df a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers' (code doesn't support missing paths, any missing paths will be represented by 0s)
#' @param paths boolean value, TRUE if you want the shortest paths returned, FALSE if you don't
#'
#' @return the number of direct patient transfers and indirect flow metrics between each facility pair. If paths = TRUE, a list of summary (pt_trans_summary) and shortest paths used (paths).
#' @export
#' @description Summarize inter-facility patient transfer network from an edge list. Direct and indirect patient flow metrics are calculated.
#' @details For more details on how patient flow is calculated, see: https://aac.asm.org/content/63/11/e01622-19.
#'
#' @examples
#' get_patient_flow(edge_df = pt_trans_df)
get_patient_flow <- function(edge_df, paths = FALSE){
  #run checks
  # could make this more general by defining column names in function, regardless of what they originally were
  check_get_patient_flow_input(edge_df = edge_df, paths = paths)

  #make pt_trans_net not factors
  edge_df$source_facil <- as.character(edge_df$source_facil)
  edge_df$dest_facil <- as.character(edge_df$dest_facil)

  #run indirect flow
  ind_flow_output <- get_indirect_flow(edge_df)
  edge_df_i <- ind_flow_output$transfer_network
  paths_list <- ind_flow_output$paths

  #remove any same-facility pairs
  pat_flow <- edge_df %>% dplyr::filter(source_facil != dest_facil)

  pt_trans_summary <- dplyr::full_join(pat_flow, edge_df_i, by = c("source_facil", "dest_facil")) %>%
    dplyr::filter(source_facil != dest_facil)

  facil_pairs <- unique(sort(sapply(1:nrow(pt_trans_summary), function(x)
    paste0(sort(c(pt_trans_summary$source_facil[x], pt_trans_summary$dest_facil[x])), collapse = ''))))

  pt_trans_summary <- lapply(facil_pairs, function(x){
    f12 <- pt_trans_summary %>% dplyr::filter(source_facil == substring(x, 1, 1) & dest_facil == substring(x, 2, 2)) %>%
      dplyr::rename(loc1 = source_facil, loc2 = dest_facil, n_transfers_f12 = n_transfers, pt_trans_metric_f12 = pt_trans_metric)
    f21 <- pt_trans_summary %>% dplyr::filter(source_facil == substring(x, 2, 2) & dest_facil == substring(x, 1, 1)) %>%
      dplyr::rename(loc1 = dest_facil, loc2 = source_facil, n_transfers_f21 = n_transfers, pt_trans_metric_f21 = pt_trans_metric)
    pt_flow_sub <- dplyr::full_join(f12, f21, by = c('loc1','loc2'))
  }) %>% dplyr::bind_rows() %>% dplyr::mutate(sum_transfers = n_transfers_f12 + n_transfers_f21,
                                sum_pt_trans_metric = pt_trans_metric_f12 + pt_trans_metric_f21)

  if(paths == TRUE){
    #return paths and summary as a list
    pt_trans_summary <- list("pt_trans_summary" = pt_trans_summary, "paths" = paths_list)
  }
  return(pt_trans_summary)
}

#' Calculate indirect patient flow from patient transfer network
#'
#' @param edge_df a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers' (code doesn't support missing paths, any missing paths will be represented by 0s)
#'
#' @return facility x facility matrix of metric of patient flow between each facility pair
#' @noRd
#'
#' @examples
#' get_indirect_flow(edge_df = pt_trans_df)
get_indirect_flow <- function(edge_df){
  #don't want to subset before getting here, need whole network for indirect
  #checks
  check_edge_df(edge_df)

  # fill in missing source and destination facilities (doesn't change results, but will error out otherwise)
  edge_df <- fill_missing_src_dest(edge_df)

  #make matrix format
  trans_mat <- tidyr::pivot_wider(edge_df, names_from = source_facil, values_from = n_transfers)
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


#' Fill in missing source and destination facilities in network edge list
#'
#' @param edge_df a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers' (code doesn't support missing paths, any missing paths will be represented by 0s)
#'
#' @return filled in edge_df
#' @noRd
#'
fill_missing_src_dest <- function(edge_df) {
  all_facils <- unique(c(as.character(edge_df$source_facil),as.character(edge_df$dest_facil)))
  not_in_source <- all_facils[!(all_facils %in% edge_df$source_facil)]
  not_in_dest <- all_facils[!(all_facils %in% edge_df$dest_facil)]
  if(length(not_in_source) != 0 | length(not_in_dest) != 0){
    edge_df$source_facil <- as.character(edge_df$source_facil)
    edge_df <- dplyr::bind_rows(edge_df,
                                dplyr::bind_cols(source_facil = not_in_source,
                                                 dest_facil = edge_df$dest_facil[1],
                                                 n_transfers = 0),
                                dplyr::bind_cols(source_facil = edge_df$source_facil[1],
                                                 dest_facil = not_in_dest,
                                                 n_transfers = 0))
  }
  edge_df <- edge_df %>% tidyr::expand(source_facil, dest_facil) %>%
    dplyr::left_join(edge_df, by = c("source_facil", "dest_facil")) %>%
    dplyr::mutate(n_transfers = ifelse(is.na(n_transfers), 0, n_transfers))
  return(edge_df)
}

