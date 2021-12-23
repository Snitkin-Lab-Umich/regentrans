#' Calculate genetic flow (Fsp)
#'
#' @inheritParams get_pair_types
#' @param fasta ape DNAbin object (i.e. from fasta file of SNPs) using read.fasta
#' @param matrix whether to output symmetric matrix (TRUE; default) or long form (FALSE)
#'
#'
#' @return facility x facility matrix with Fsp values
#' @export
#' @details Genetic flow (Fsp) is described in Donker et al. 2017
#' (mgen.microbiologyresearch.org/pubmed/content/journal/mgen/10.1099/mgen.0.000113).
#' Only bi-allelic sites are included when computing Fsp.
#' The Fsp values are between 0 and 1 where lower values indicate more similar populations.
#' Note that the current implementation of this function is fairly slow, visit https://github.com/nateosher/RPTfast for a faster implementation
#'
#' @examples
#' \dontrun{
#' # This takes a long time to run right now!
#' locs <- metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' pt <- metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()
#' facil_fsp <- get_genetic_flow(aln, locs, matrix = TRUE, pt)
#' }


#######################try to make sapply########
get_genetic_flow <- function(fasta, locs, matrix = TRUE, pt){
  #check the DNAbin object and locs
  check_facility_fsp_input(fasta, locs, matrix, pt)
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
  #unique list of locations
  locs_unique <- unique(unname(locs_subset))

  #generate metasequences if pt info is given
  if(!is.null(pt)){
    fasta_sub_meta <- make_meta_seqs(fasta_sub, locs_subset, pt[isolates])
  }
  else{
    fasta_sub_meta <- fasta_sub
  }

  #TO DO:re-subset the locs df based on the now rownames of the fasta sub meta
  locs_subset <- locs_subset[rownames(fasta_sub_meta)]
  #order the locs object to match the order of rownames of the fasta
  #this might be redundant
  locs_subset <- locs_subset[order(match(names(locs_subset),rownames(fasta_sub_meta)))]

  #change the rownames of the fasta to the location names
  rownames(fasta_sub_meta) <- unname(locs_subset)
  #make a list of the location names
  sample_locs <- rownames(fasta_sub_meta)


  #CALCULATE INTRA- AND INTER-FACILITY DISTANCE
  facil_dist <- data.frame(sapply(locs_unique, function(f1){
    sapply(locs_unique, function(f2){
      if (f1 == f2) {return(0)}
      #subset fasta file to just that those locations
      subset_snp_mat = fasta_sub_meta[sample_locs %in% c(f1, f2), ]
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
#' @inheritParams get_genetic_flow
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
#' @inheritParams get_genetic_flow
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

make_meta_seqs <- function(fasta, locs, pt){
  #TO DO: add checks
  #remake metadata df with location, pt, and isolate ID
  combos <- as.data.frame(cbind(pt, names(pt))) %>% dplyr::left_join(as.data.frame(cbind(locs, names(locs))))
  #subset to single pt/locs combos
  combos_2 <- combos %>% dplyr::distinct(pt, locs, .keep_all = TRUE)
  #if there are not multiple patients from the same location, just return the normal fasta
  if(nrow(combos) == nrow(combos_2)){
    return(fasta)
  }
  #otherwise
  #Find major allele at each position, make a “reference”
  ref <- find_major_alleles(as.character(fasta))
  #make df to count number of unique combos
  combos_3 <- combos %>% dplyr::count(pt, locs)
  #make fasta character
  fasta_char <- as.character(fasta)
  #across all patient location combos
  fasta_subs <- data.frame(t(apply(combos_3,1,function(x){
    #find the isolate ID
    ID <- combos %>% dplyr::filter(pt == x[1][[1]], locs == x[2][[1]]) %>% dplyr::select(V2)
    #if there is only one unique pair, return that sequence as is
    if(x[3] == 1){
      metasequence <- fasta_char[rownames(fasta_char) == ID[[1]], ]
    }
    #if there are more than one, make a metasequence
    else{
      metasequence <- find_major_alleles(fasta_char[rownames(fasta_char) %in% ID[[1]], ], ref = ref)
    }
    #return one ID so we can map it back to the location and the sequence
    return(c(ID[[1]][1], metasequence))

  })))
  #make first col (seq_ID) colnames
  rownames(fasta_subs) <- fasta_subs[,1]
  #remove that column
  fasta_subs <- fasta_subs[,-1]
  #return the metasequence as a DNAbin
  return(ape::as.DNAbin(as.matrix(fasta_subs)))
}

find_major_alleles <- function(fasta, ref = NULL){
  #make character fasta (must be entered as character)
  #TO DO: write tests for fasta and ref as character
  fasta_2 <- fasta
  #if there isn't a ref provided, we are making the ref
  if(is.null(ref)){
    #return the major allele (most common) at each position
    ref <- apply(fasta_2, 2, FUN = function(x){
      #only problem is if there is a tie, which is less of a problem with larger sample size
      names(which.max(table(x)))
    })
  }
  #if there is already a ref, we are finding the consensus sequence for a patient
  else{
    #apply across all positions for the genome
    ref <- sapply(1:ncol(fasta_2), FUN = function(i){
      #subset to that column (position)
      x = fasta_2[,i]
      #make a table of frequency of alleles in that seq
      tab = table(x)
      #check if there are any values differing from the reference allele
      #if there is only one allele, return it
      if(length(tab) == 1){
        return(names(tab))
      }
      #if there is more than one value, return the one that isn't the reference allele
      else{
        #find other values
        alt_allele <- names(tab)[names(tab) != ref[i]]
        #return first other value
        return(alt_allele[1])
      }
    })

  }

  return(ref)
}

