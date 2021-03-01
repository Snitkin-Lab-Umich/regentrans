#Joyce's code (joycewang914 on github)
# from joyce. calculates fst as described in Donker et al. 2017
# input:
#   - snp_dist = ape DNAbin object (i.e. from fasta file of SNPs) using read.fasta?
#   - facil = facility names of interest
# returns fst for each pair of facilities


#' Make facility x facility matrix with Fsp values. Calculates fst as described in Donker et al. 2017
#'
#' @param snp_dist ape DNAbin object (i.e. from fasta file of SNPs) using read.fasta
#' @param locs locations names for pairwise comparison
#'
#' @return matrix of facility x facility matrix with Fsp values. Only bi-allelic sites. Fsp values bween 0 (HP=HS) and 1 (Hp = 0)
#' @export
#'
#' @examples



#######################try to make sapply########
get_facility_fsp <- function(fasta, locs){
  #check the DNAbin object and locs
  check_facility_fsp(fasta, locs)
  #TODO: add check to make sure at least 2 of the locs appear more than once
  #make a vector of only locs that appear more than once
  locs_over_one <- which(unlist(table(locs) > 1))
  locs_subset <- locs[locs %in% names(locs_over_one)]
  #make a list of the ones they have in common for subsetting
  isolates <- intersect(names(locs_subset), rownames(fasta))
  #subset the DNAbin object to the samples they have in common
  fasta<-fasta[isolates,]
  #subset the locs object to the samples they have in common
  locs_subset <- locs_subset[isolates]
  #order the locs object to match the order of rownames of the fasta
  locs_subset <- locs_subset[order(match(names(locs_subset),rownames(fasta)))]
  #change the rownames of the fasta to the location names
  rownames(fasta) <- unname(locs_subset)
  #make a list of the location names
  sample_locs = rownames(fasta)
  #unique list of locations
  locs_unique = unique(unname(locs_subset))

  #CALCULATE INTRA- AND INTER-FACILITY DISTANCE
  facil_dist <- data.frame(sapply(locs_unique, function(f1){
    sapply(locs_unique, function(f2){
      if (f1 == f2) {return(0)}
      #subset fasta file to just that those locations
      subset_snp_mat = fasta[sample_locs %in% c(f1, f2), ]
      #make sure the position is not all unknown or no variance, subset to the ones that have some variation
      subset_snp_mat = subset_snp_mat[,apply(subset_snp_mat, 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
      #figure out which are from each facility
      subset_f1 = rownames(subset_snp_mat) %in% f1
      subset_f2 = rownames(subset_snp_mat) %in% f2
      #BETWEEN POPLUATION VARIATION
      #loop over columns (variants)
      between = apply(subset_snp_mat, 2, FUN = function(x){
        #get alleles present at the site
        alleles = names(table(as.character(x)))
        #skip multi-allelic sites
        if (length(alleles) > 2){0} else{
          #find allele frequency for each allele at each site
          #make function to call x4
          f1_allele1 = sum(as.character(x)[subset_f1] %in% alleles[1])/sum(subset_f1)
          f1_allele2 = sum(as.character(x)[subset_f1] %in% alleles[2])/sum(subset_f1)
          f2_allele1 = sum(as.character(x)[subset_f2] %in% alleles[1])/sum(subset_f2)
          f2_allele2 = sum(as.character(x)[subset_f2] %in% alleles[2])/sum(subset_f2)
          f1_allele1 * f1_allele2 * f2_allele1 * f2_allele2}
      })
      #sum
      between_sum = sum(between)
      #WITHIN POPULATION 1 VARIATION
      #subset to first location
      f1_subset_snp_mat = subset_snp_mat[subset_f1,apply(subset_snp_mat[subset_f1,], 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
      #same as above
      within_f1 = apply(f1_subset_snp_mat, 2, FUN = function(x){
        alleles = names(table(as.character(x)))
        if (length(alleles) > 2){0}else{
          f1_allele1 = sum(as.character(x) %in% alleles[1])/length(x)
          f1_allele2 = sum(as.character(x) %in% alleles[2])/length(x)
          (f1_allele1 * f1_allele2)^2
        }
      })
      within_f1_sum = sum(within_f1)
      #WITHIN POPULATION 2 VARIATION
      f2_subset_snp_mat = subset_snp_mat[subset_f2,apply(subset_snp_mat[subset_f2,], 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
      within_f2 = apply(f2_subset_snp_mat, 2, FUN = function(x){
        alleles = names(table(as.character(x)))
        if (length(alleles) > 2){0}else{
          f2_allele1 = sum(as.character(x) %in% alleles[1])/length(x)
          f2_allele2 = sum(as.character(x) %in% alleles[2])/length(x)
          (f2_allele1 * f2_allele2)^2}
      })
      within_f2_sum = sum(within_f2)
      Fsp = (((within_f1_sum + within_f2_sum) / 2) - between_sum) / ((within_f1_sum + within_f2_sum) / 2)
      return(Fsp)
    })#end loop 1
  }))#end loop 2
  #add row and column names
  rownames(facil_dist) <- locs_unique
  colnames(facil_dist) <- locs_unique
  return(facil_dist);
}#end facility_fst


########################3normal!#############
# get_facility_fsp <- function(fasta, locs){
#   #check the DNAbin object and locs
#   check_facility_fsp(fasta, locs)
#   #make a list of the ones they have in common for subsetting
#   isolates <- intersect(names(locs), rownames(fasta))
#   #subset the DNAbin object to the samples they have in common
#   fasta<-fasta[isolates,]
#   #subset the locs object to the samples they have in common
#   locs <- locs[isolates]
#   #order the locs object to match the order of rownames of the fasta
#   locs <- locs[order(match(names(locs),rownames(fasta)))]
#   #change the rownames of the fasta to the location names
#   rownames(fasta) <- unname(locs)
#   #make a list of the location names
#   sample_locs = rownames(fasta)
#   #unique list of locations
#   locs_unique = unique(unname(locs))
#
#   #CALCULATE INTRA- AND INTER-FACILITY DISTANCE
#   facil_dist = matrix(0, ncol = length(locs_unique), nrow = length(locs_unique), dimnames = list(locs_unique, locs_unique))
#   for(f1 in locs_unique){
#     for(f2 in locs_unique){
#       if (f1 == f2) {next}
#       #subset fasta file to just that those locations
#       subset_snp_mat = fasta[sample_locs %in% c(f1, f2), ]
#       #make sure the position is not all unknown or no variance, subset to the ones that have some variation
#       subset_snp_mat = subset_snp_mat[,apply(subset_snp_mat, 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
#       #figure out which are from each facility
#       subset_f1 = rownames(subset_snp_mat) %in% f1
#       subset_f2 = rownames(subset_snp_mat) %in% f2
#       #BETWEEN POPLUATION VARIATION
#       #loop over columns (variants)
#       between = apply(subset_snp_mat, 2, FUN = function(x){
#         #get alleles present at the site
#         alleles = names(table(as.character(x)))
#         #skip multi-allelic sites
#         if (length(alleles) > 2){0} else{
#           #find allele frequency for each allele at each site
#           #make function to call x4
#           f1_allele1 = sum(as.character(x)[subset_f1] %in% alleles[1])/sum(subset_f1)
#           f1_allele2 = sum(as.character(x)[subset_f1] %in% alleles[2])/sum(subset_f1)
#           f2_allele1 = sum(as.character(x)[subset_f2] %in% alleles[1])/sum(subset_f2)
#           f2_allele2 = sum(as.character(x)[subset_f2] %in% alleles[2])/sum(subset_f2)
#           f1_allele1 * f1_allele2 * f2_allele1 * f2_allele2}
#       })
#       #sum
#       between_sum = sum(between)
#       #WITHIN POPULATION 1 VARIATION
#       #subset to first location
#       f1_subset_snp_mat = subset_snp_mat[subset_f1,apply(subset_snp_mat[subset_f1,], 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
#       #same as above
#       within_f1 = apply(f1_subset_snp_mat, 2, FUN = function(x){
#         alleles = names(table(as.character(x)))
#         if (length(alleles) > 2){0}else{
#           f1_allele1 = sum(as.character(x) %in% alleles[1])/length(x)
#           f1_allele2 = sum(as.character(x) %in% alleles[2])/length(x)
#           (f1_allele1 * f1_allele2)^2
#         }
#       })
#       within_f1_sum = sum(within_f1)
#       #WITHIN POPULATION 2 VARIATION
#       f2_subset_snp_mat = subset_snp_mat[subset_f2,apply(subset_snp_mat[subset_f2,], 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
#       within_f2 = apply(f2_subset_snp_mat, 2, FUN = function(x){
#         alleles = names(table(as.character(x)))
#         if (length(alleles) > 2){0}else{
#           f2_allele1 = sum(as.character(x) %in% alleles[1])/length(x)
#           f2_allele2 = sum(as.character(x) %in% alleles[2])/length(x)
#           (f2_allele1 * f2_allele2)^2}
#       })
#       within_f2_sum = sum(within_f2)
#       Fsp = (((within_f1_sum + within_f2_sum) / 2) - between_sum) / ((within_f1_sum + within_f2_sum) / 2)
#       facil_dist[f1,f2] = Fsp
#     }#end loop 1
#   }#end loop 2
#     return(facil_dist);
# }#end facility_fst

############################code from joyce's github###############################################################################
  # facility_fsp -
  # Input:
  #   1) snp_dist = variant matrix (sample genomes x single nucleotide variants)
  #   2) facil = facilties for pairwise comparison
  #      Note that our sample genomes are named as "3425-5015-0-R2X" where:
  #     "5015" = patient ID
  #     "5" (first character of patient ID) = facility ID
  # Returns a facil x facil matrix with Fsp values

  # get_facility_fsp <- function(snp_dist, facil)
  # {
  #
  #   # SUBSET SAMPLES BASED ON FACILITY
  #   #gets list of facility IDs from the rownames of the snp_dist matrix
  #   sample_facil = substr(sapply(rownames(snp_dist), FUN = function(x){strsplit(x, "-")[[1]][2]}), 1, 1)
  #
  #   #CALCULATE INTRA- AND INTER-FACILITY DISTANCE
  #   #makes empty data frame of facility x facility
  #   facil_dist = matrix(0, ncol = length(facil), nrow = length(facil), dimnames = list(facil, facil))
  #
  #   for(f1 in facil)
  #   {
  #
  #     for(f2 in facil)
  #     {
  #
  #       if (f1 == f2) {next}
  #
  #       subset_snp_mat = snp_dist[sample_facil %in% c(f1, f2), ]
  #       subset_snp_mat = subset_snp_mat[,apply(subset_snp_mat, 2, FUN = function(x){sum(x != x[1] | x == 'n') > 0})]
  #       subset_f1 = substr(sapply(rownames(subset_snp_mat), FUN = function(x){strsplit(x, "-")[[1]][2]}), 1, 1) %in% f1
  #       subset_f2 = substr(sapply(rownames(subset_snp_mat), FUN = function(x){strsplit(x, "-")[[1]][2]}), 1, 1) %in% f2
  #
  #       #BETWEEN POPLUATION VARIATION
  #
  #       between = apply(subset_snp_mat, 2, FUN = function(x){
  #
  #         alleles = names(table(as.character(x)))
  #
  #         if (length(alleles) > 2){0} else{
  #
  #           f1_allele1 = sum(as.character(x)[subset_f1] %in% alleles[1])/sum(subset_f1)
  #           f1_allele2 = sum(as.character(x)[subset_f1] %in% alleles[2])/sum(subset_f1)
  #
  #           f2_allele1 = sum(as.character(x)[subset_f2] %in% alleles[1])/sum(subset_f2)
  #           f2_allele2 = sum(as.character(x)[subset_f2] %in% alleles[2])/sum(subset_f2)
  #
  #           f1_allele1 * f1_allele2 * f2_allele1 * f2_allele2
  #         }
  #       })
  #
  #       between_sum = sum(between)
  #
  #       #WITHIN POPULATION 1 VARIATION
  #
  #       f1_subset_snp_mat = subset_snp_mat[subset_f1,apply(subset_snp_mat[subset_f1,], 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
  #
  #       within_f1 = apply(f1_subset_snp_mat, 2, FUN = function(x){
  #
  #         alleles = names(table(as.character(x)))
  #
  #         if (length(alleles) > 2){0}else{
  #
  #           f1_allele1 = sum(as.character(x) %in% alleles[1])/length(x)
  #           f1_allele2 = sum(as.character(x) %in% alleles[2])/length(x)
  #
  #           (f1_allele1 * f1_allele2)^2
  #         }
  #       })
  #
  #       within_f1_sum = sum(within_f1)
  #
  #
  #       #WITHIN POPULATION 2 VARIATION
  #
  #       f2_subset_snp_mat = subset_snp_mat[subset_f2,apply(subset_snp_mat[subset_f2,], 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]
  #
  #       within_f2 = apply(f2_subset_snp_mat, 2, FUN = function(x){
  #
  #         alleles = names(table(as.character(x)))
  #
  #         if (length(alleles) > 2){0}else{
  #
  #           f2_allele1 = sum(as.character(x) %in% alleles[1])/length(x)
  #           f2_allele2 = sum(as.character(x) %in% alleles[2])/length(x)
  #
  #           (f2_allele1 * f2_allele2)^2}
  #       })
  #
  #       within_f2_sum = sum(within_f2)
  #
  #       Fsp = (((within_f1_sum + within_f2_sum) / 2) - between_sum) / ((within_f1_sum + within_f2_sum) / 2)
  #
  #
  #       facil_dist[f1,f2] = Fsp
  #
  #     }#end for
  #
  #
  #   }#end for
  #
  #   return(facil_dist);
  #
  # }#end facility_fsp
  #
  #

#### notes##########
# x <- 1:5
# y <- 1:5
# for(i in x){
#   for(j in y){
#     print(i+j)
#   }
# }
#
# #makes the matrix
# sapply(x, function(i){
#   sapply(y, function(j){
#     if(i == j) {return()}
#     else{return(i+j)}
#     i+j
#   })
# })
#
# test <- data.frame(sapply(x, function(i){
#   sapply(y, function(j){
#     if(i == j) {return()}
#     else{return(i+j)}
#     i+j
#   })
# }))
#
#
# #make small fasta with a few to check
# locs <- locs[1:100]
