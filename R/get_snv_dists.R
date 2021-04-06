#' Get matrix of inter- and intra-facility pariwise snv distances
#'
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patient that isolate was taken from with the name being sample ID (optional)
#' @param pt_trans_net a dataframe representing a patient transfer network of 3 cols: 'source_facil', 'dest_facil, and 'n_transfers'
#'
#' @return snp_facility_paris, a matrix of inter- and intra-facility pariwise snv distances
#' @export
#'
#' @examples to reduce this dataframe to include only one row that represents each pair, use the subset_pairs() function on the output of this function
get_snv_dists <- function(dists, locs, pt = NULL, pt_trans_net = NULL){
  #checks
  check_get_snv_dists_input(dists, locs, pt, pt_trans_net)

  #check pt if it isn't missing
  if(!is.null(pt)){
    #subset pt and locs to represent the same isolates,
    #and make sure we only subset to the rownames dists has in common
    isolates <- intersect(intersect(names(locs), names(pt)), rownames(dists))
    #subset pt
    pt_sub <- pt[isolates]
  }
  else{
    #make the subsetted isolates object if there is no pt
    isolates <- intersect(names(locs), rownames(dists))
  }

  #subset by locs
  #list ones in common before subsetting
  loc_sub <- locs[isolates]

  #subset dists to isolates
  dists_sub <- dists[isolates, isolates]

  #make df
  snps <- na.omit(data.frame(as.table(as.matrix(dists_sub))))
  #change freq colname?
  colnames(snps) <- c("Isolate1", "Isolate2", "Pairwise_Dists")

  #add locs
  snps$Loc1 <- loc_sub[snps$Isolate1]
  snps$Loc2 <- loc_sub[snps$Isolate2]

  #do pt stuff
  if(!is.null(pt)){
    #add pts
    snps$Patient1 <- pt_sub[snps$Isolate1]
    snps$Patient2 <- pt_sub[snps$Isolate2]
    #add labels
    snp_facility_pairs <- dplyr::bind_cols(snps %>% filter(Patient1 != Patient2) %>% mutate(Pair_Type=ifelse(Loc1==Loc2,'Intra-facility pair','Inter-facility pair')))
  }
  else{
    #add labels
    snp_facility_pairs <- dplyr::bind_cols(snps %>% filter(Isolate1 != Isolate2) %>% mutate(Pair_Type=ifelse(Loc1==Loc2,'Intra-facility pair','Inter-facility pair')))
  }
  #if there is a patient transfer network, add it
  if(!is.null(pt_trans_net)){
    #edit the network to make sure pairs exist forward and backward
    #make a dataframe?
    tidy_pt_trans_net <- tidyr::pivot_wider(pt_trans_net, names_from = source_facil, values_from = n_transfers)
    #if the matrix is symmetric
    #symmetry <- isSymmetric(as.matrix(tidy_pt_trans_net[,2:ncol(tidy_pt_trans_net)]))
    pt_trans_net <- data.frame(as.table(as.matrix(tidy_pt_trans_net))) %>% filter(Var2 != "dest_facil") %>% `colnames<-`(c("source_facil", "dest_facil", "n_transfers"))
    snp_facility_pairs <- snp_facility_pairs %>% dplyr::left_join(pt_trans_net, by = c("Loc1" = "source_facil", "Loc2" = "dest_facil"))
    snp_facility_pairs$n_transfers <- as.numeric(snp_facility_pairs$n_transfers)
    #edit the ones that didn't match to be something
  }
  #return snp matrix
  return(snp_facility_pairs)
}
