#' Calculate gene flow (Fsp)
#'
#' @inheritParams get_snv_dists
#' @param fasta ape DNAbin object (i.e. from fasta file of SNPs) using read.fasta
#' @param matrix whether to output symmetric matrix (TRUE; default) or long form (FALSE)
#'
#' @return facility x facility matrix with Fsp values
#' @export
#' @details Fsp is described in Donker et al. 2017
#' (mgen.microbiologyresearch.org/pubmed/content/journal/mgen/10.1099/mgen.0.000113).
#' Only bi-allelic sites are included when computing Fsp.
#' The Fsp values are between 0 and 1 where lower values indicate more similar populations.
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
