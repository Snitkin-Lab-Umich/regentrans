#' Make facility x facility matrix with fsp values. Calculates fst as described in Donker et al. 2017
#'
#' @param fasta ape DNAbin object (i.e. from fasta file of SNPs) using read.fasta
#' @param locs locations names for pairwise comparison
#' @param matrix whether to output symmetric matrix (TRUE; default) or long form (FALSE)
#'
#' @return matrix of facility x facility matrix with fsp values. Only bi-allelic sites. fsp values bween 0 (HP=HS) and 1 (Hp = 0)
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' facil_fsp <- get_facility_fsp(aln, locs, matrix = TRUE)
#' }


#######################try to make sapply########
get_facility_fsp <- function(fasta, locs, matrix = TRUE){
  #check the DNAbin object and locs
  check_facility_fsp_input(fasta, locs, matrix)
  #make a vector of only locs that appear more than once
  locs_over_one <- which(unlist(table(locs) > 1))
  locs_subset <- locs[locs %in% names(locs_over_one)]
  #make a list of the ones they have in common for subsetting
  isolates <- intersect(names(locs_subset), rownames(fasta))
  #subset the DNAbin object to the samples they have in common
  fasta_sub<-fasta[isolates,]
  #subset the locs object to the samples they have in common
  locs_subset <- locs_subset[isolates]
  #order the locs object to match the order of rownames of the fasta
  #this might be redundant
  locs_subset <- locs_subset[order(match(names(locs_subset),rownames(fasta_sub)))]
  #change the rownames of the fasta to the location names
  rownames(fasta_sub) <- unname(locs_subset)
  #make a list of the location names
  sample_locs <- rownames(fasta_sub)
  #unique list of locations
  locs_unique <- unique(unname(locs_subset))
  #CALCULATE INTRA- AND INTER-FACILITY DISTANCE
  facil_dist <- data.frame(sapply(locs_unique, function(f1){
    sapply(locs_unique, function(f2){
      if (f1 == f2) {return(0)}
      #subset fasta file to just that those locations
      subset_snp_mat = fasta_sub[sample_locs %in% c(f1, f2), ]
      #make sure the position is not all unknown or no variance, subset to the ones that have some variation
      subset_snp_mat = subset_snp_mat[,apply(subset_snp_mat, 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
      #figure out which are from each facility
      subset_f1 = rownames(subset_snp_mat) %in% f1
      subset_f2 = rownames(subset_snp_mat) %in% f2
      #BETWEEN POPLUATION VARIATION
      #for each position
      between = apply(subset_snp_mat, 2, FUN = function(x){
        #get alleles present at the site
        alleles = names(table(as.character(x)))
        #skip multi-allelic sites
        if (length(alleles) > 2){0} else{
          #find allele frequency for each allele at each site
          f1_allele1 = get_allele_freq_btwn(x = x, subset = subset_f1, allele_n = 1, alleles = alleles)
          f1_allele2 = get_allele_freq_btwn(x, subset_f1, 2, alleles)
          f2_allele1 = get_allele_freq_btwn(x, subset_f2, 1, alleles)
          f2_allele2 = get_allele_freq_btwn(x, subset_f2, 2, alleles)
          #calculate between pop variation for each allele site?
          f1_allele1 * f1_allele2 * f2_allele1 * f2_allele2}
      })
      #sum
      between_sum = sum(between)
      #WITHIN POPULATION 1 VARIATION
      within_f1_sum <- get_within_pop_var(subset_snp_mat, subset_f1)
      #WITHIN POPULATION 2 VARIATION
      within_f2_sum <- get_within_pop_var(subset_snp_mat, subset_f2)
      #calculate fsp
      fsp = (((within_f1_sum + within_f2_sum) / 2) - between_sum) / ((within_f1_sum + within_f2_sum) / 2)
      return(fsp)
    })#end loop 1
  }))#end loop 2
  #add row and column names
  rownames(facil_dist) <- locs_unique
  colnames(facil_dist) <- locs_unique
  #change to long form if that is specified
  if(!matrix){ facil_dist <- make_long_form(facil_dist) }
  return(facil_dist);
}#end facility_fst


#' Find allele frequency for each allele at each site
#'
#' @inheritParams get_facility_fsp
#'
#' @noRd
#'
get_allele_freq_btwn <- function(x, subset, allele_n, alleles){
  #checks
  check_allele_freq_input(x, subset, allele_n, alleles)
  #calculate
  return(sum(as.character(x)[subset] %in% alleles[allele_n])/sum(subset))
}

get_allele_freq_within <- function(x, allele_n, alleles){
  #checks
  check_allele_freq_input(x, subset = NULL, allele_n, alleles)
  #calculate
  return(sum(as.character(x) %in% alleles[allele_n])/length(x))
}

#' Get within population variation
#'
#' @inheritParams get_facility_fsp
#'
#' @noRd
#'
get_within_pop_var <- function(subset_snp_mat, subset){
  #checks
  check_within_pop_var_inputs(subset_snp_mat, subset)

  f_subset_snp_mat = subset_snp_mat[subset,apply(subset_snp_mat[subset,], 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
  within_f = apply(f_subset_snp_mat, 2, FUN = function(x){
    alleles = names(table(as.character(x)))
    if (length(alleles) > 2){0}else{
      f_allele1 = get_allele_freq_within(x, 1, alleles)
      f_allele2 = get_allele_freq_within(x, 2, alleles)

      (f_allele1 * f_allele2)^2}
  })
  return(sum(within_f))
}

#' Make symmetric matrix long-form where each row represents a cell in the matrix
#'
#' @param facil_dist symmetric matrix that you want to convert to long form
#'
#' @return long form data where each row represents a cell in the matrix (rowname, column name, value)
#' @export
#'
#' @examples
#' make_long_form(fsp)
make_long_form <- function(facil_dist, col_names = c('loc1', 'loc2', 'fsp')){
  #check that it is a symmetric matrix
  check_long_form_input(facil_dist, col_names)
  #change to longform
  facil_dist_long <- stats::na.omit(data.frame(as.table(as.matrix(facil_dist)))) %>% dplyr::filter(Freq != 0)
  colnames(facil_dist_long) <- col_names

  facil_pairs <- sapply(1:nrow(facil_dist_long), function(x)
    paste0(sort(c(as.character(facil_dist_long$loc1[x]), as.character(facil_dist_long$loc2[x]))), collapse = ''))

  facil_dist_long$loc1 <- sapply(facil_pairs, function(x) substring(x, 1, 1))
  facil_dist_long$loc2 <- sapply(facil_pairs, function(x) substring(x, 2, 2))

  return(subset_pairs(facil_dist_long))
}


