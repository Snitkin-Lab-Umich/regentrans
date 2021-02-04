#Joyce's code (joycewang914 on github)

# facility_fsp -
# Input:
#   1) snp_dist = variant matrix (sample genomes x single nucleotide variants)
#   2) facil = facilties for pairwise comparison
#      Note that our sample genomes are named as "3425-5015-0-R2X" where:
#     "5015" = patient ID
#     "5" (first character of patient ID) = facility ID
# Returns a facil x facil matrix with Fsp values

facility_fsp <- function(snp_dist, facil)
{

  # SUBSET SAMPLES BASED ON FACILITY
  sample_facil = substr(sapply(rownames(snp_dist), FUN = function(x){strsplit(x, "-")[[1]][2]}), 1, 1)

  #CALCULATE INTRA- AND INTER-FACILITY DISTANCE
  facil_dist = matrix(0, ncol = length(facil), nrow = length(facil), dimnames = list(facil, facil))

  for(f1 in facil)
  {

    for(f2 in facil)
    {

      if (f1 == f2) {next}

      subset_snp_mat = snp_dist[sample_facil %in% c(f1, f2), ]
      subset_snp_mat = subset_snp_mat[,apply(subset_snp_mat, 2, FUN = function(x){sum(x != x[1] | x == 'n') > 0})]
      subset_f1 = substr(sapply(rownames(subset_snp_mat), FUN = function(x){strsplit(x, "-")[[1]][2]}), 1, 1) %in% f1
      subset_f2 = substr(sapply(rownames(subset_snp_mat), FUN = function(x){strsplit(x, "-")[[1]][2]}), 1, 1) %in% f2

      #BETWEEN POPLUATION VARIATION

      between = apply(subset_snp_mat, 2, FUN = function(x){

        alleles = names(table(as.character(x)))

        if (length(alleles) > 2){0} else{

          f1_allele1 = sum(as.character(x)[subset_f1] %in% alleles[1])/sum(subset_f1)
          f1_allele2 = sum(as.character(x)[subset_f1] %in% alleles[2])/sum(subset_f1)

          f2_allele1 = sum(as.character(x)[subset_f2] %in% alleles[1])/sum(subset_f2)
          f2_allele2 = sum(as.character(x)[subset_f2] %in% alleles[2])/sum(subset_f2)

          f1_allele1 * f1_allele2 * f2_allele1 * f2_allele2
        }
      })

      between_sum = sum(between)

      #WITHIN POPULATION 1 VARIATION

      f1_subset_snp_mat = subset_snp_mat[subset_f1,apply(subset_snp_mat[subset_f1,], 2, FUN = function(x){sum(x != x[1] | x == 'N') > 0})]

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


      facil_dist[f1,f2] = Fsp

    }#end for


  }#end for

  return(facil_dist);

}#end facility_fsp
