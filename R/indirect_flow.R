#' Calculate indirect patient flow from patient transfer network
#'
#' @param pt_trans_net a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers'
#'
#' @return facility x facility matrix of metric of patient flow between each facility pair
#' @export
#'
#' @examples indirect_flow(pt_trans_net)
indirect_flow <- function(pt_trans_net){
  #don't want to subset before getting here, need whole network for indirect
  #checks
  check_pt_trans_net(pt_trans_net, unique(c(pt_trans_net$source_facil, pt_trans_net$dest_facil)))

  #make matrix format
  trans_mat <- tidyr::pivot_wider(pt_trans_net, names_from = source_facil, values_from = n_transfers)
  trans_mat <- as.data.frame(trans_mat[,2:ncol(trans_mat)])
  rownames(trans_mat) <- colnames(trans_mat)

  #make graph
  g <- graph_from_adjacency_matrix(as.matrix(trans_mat),mode='directed',weighted = TRUE)


  #do some magic??

  #name nodes in network that we have data for

  #find shortest path function -> igraph::shortest.paths()

  #modify edge weights from n facilities -> normalize/invert

  #find shortest path between all nodes of interest

  #square matrix of facilities of interest
  #metric of patient flow between them


  #will this include direct, overwrite the other one or have to be added to the other one?

}

##########################code from zena########################
# library(igraph)
# library(readxl)
# source('../../../Project_Penn_KPC/Analysis/Patient_transfer/lib/pt_transf_network_functions_v2.R')
#
# ptmat <- read_excel('../2020-11-10_facilities/data/Matrices LOS 09292016.xlsx',sheet = 1)
# r_fids <- structure(c("N_vsnf", "C_ltach", "D_ltach", "G_icu", "D_icu",
#                       "Y_icu", "J_vsnf", "Q_icu", "H_vsnf", "M_vsnf", "S_icu",
#                       "F_ltach", "L_vsnf", "B_ltach", "A_ltach", "G_ltach", "I_vsnf",
#                       "K_vsnf", "O_vsnf", "E_ltach", "M_icu", "L_icu"),
#                     .Names = c("ALDE_820_V", "THC_2544_L", "THC_6130_L", "MERC_2525_H",
#                                "LORE_645_H", "COMM_5645_H", "GLEN_22660_V", "NORT_251_H",
#                                "GLEN_2451_V", "STA_1725_V", "THE_1740_H", "ADVO_3435_L", "ELMW_7733_V",
#                                "THC_365_L", "THC_4058_L", "RML_5601_L", "GLEN_8333_V", "BALL_9300_V",
#                                "OAK_9525_V", "PRES_100_L",  "SAIN_2875_H", "NORW_1044_H"))
# names(ptmat) <- sapply(names(ptmat), function(x){
#   if(x %in% paste0('To_',names(r_fids))){
#     nam <- paste0('To_',r_fids[gsub('To_','',x)])
#   }else{
#     x
#   }
# })
#
# ptmat$UNIQUE_ID <- sapply(ptmat$UNIQUE_ID, function(x){
#   if(x %in% names(r_fids)){
#     nam <- r_fids[x]
#   }else{
#     x
#   }
# })
#
# sub <- data.frame(ptmat[ptmat$UNIQUE_ID %in% r_fids,names(ptmat) %in% c('UNIQUE_ID',paste0('To_',r_fids))])
# rownames(sub) <- sub$UNIQUE_ID
#
# rn_ptmat <- ptmat$UNIQUE_ID
# ptmat <- as.matrix(ptmat[,2:ncol(ptmat)])
# rownames(ptmat) <- rn_ptmat
#
# g <- graph_from_adjacency_matrix(as.matrix(ptmat),mode='directed',weighted = TRUE)
#
# V(g)$size = degree(g, mode="all")*0.01
# E(g)$width = E(g)$weight/500
#
# V(g)$labels = sapply(1:length(V(g)$size),function(x) {
#   if(!names(V(g))[x] %in% paste0('To_',r_fids)){
#     NA
#   }else{
#     (r_fids)[paste0('To_',r_fids) == names(V(g))[x]]
#   }
# })
#
# V(g)$facil_type = unlist(sapply(V(g)$labels, function(x){
#   if(is.na(x)){
#     NA
#   }else{
#     strsplit(x,split='_')[[1]][2]
#   }
# }))
#
# # hist(E(simplify(g))$weight,breaks=100)
# # hist(E(simplify(g))$weight[E(simplify(g))$weight < 10],breaks=100)
# # hist(E(simplify(g))$weight[E(simplify(g))$weight >= 2],breaks=100)
# # g = delete_edges(g, E(g)[E(g)$weight<=10])
#
# grps = components(g, mode='weak')
# g = delete_vertices(g,names(grps$membership)[grps$membership!=1])
#
# # plot(simplify(g),vertex.label=V(g)$labels,edge.arrow.size=0.05,layout=layout_(g,with_fr()),
# #          vertex.label.cex=0.5,vertex.frame.color=NA,vertex.color='orange',
# #          edge.color='grey',edge.width=E(g)$width/1000)
#
# out_strength = strength(g,mode='out') # get number of outgoing patient transfers for each vertex
# tail_vert = tail_of(g,E(g)) # get tail (source) vertex for each edge
# edwt_sum = sapply(names(tail_vert), function(x) out_strength[names(out_strength) == x]) # get number of outgoing patient transfers of tail vertex for each edge
# E(g)$weight_norm = -log10(E(g)$weight/edwt_sum) # normalize edge weight by number of outgoing patient transfers of source vertex and take negative log (to use to calculate shortest paths)
# E(g)$weight_norm_nolog = E(g)$weight/edwt_sum
#
# # GET VERTEX NAMES (UNIQUE NUMBERS) AND ALIASES FOR LTACHS OF INTEREST
# ltach_vertex_inds = as.vector(V(g)[!is.na(V(g)$labels)])
# names(ltach_vertex_inds) = V(g)$labels[!is.na(V(g)$labels)]
#
# g = set.vertex.attribute(g, "name", value=V(g))
#
# save(g, file = paste0(format(Sys.time(), "%Y-%m-%d"),'_graph.RData'))
#
# num_k = 100
# ksp = get_dist_mat_ksp(num_k, g, ltach_vertex_inds, weight = E(g)$weight_norm)
# save(ksp,file = paste0(format(Sys.time(), "%Y-%m-%d"),'_ksp',num_k,'.RData'))


