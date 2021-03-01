# some notes:
# we might want to change the function names (should start with a verb, be informative and not too long)
# we might want to change the format/colnames of the output as well
# good to include checks to see if the input is correct (within the function; can have a separate script for all check functions) and unit tests to check the functions (also in a separate script; we'll hopefully find bugs this way!)
# we'll want to make some arguments optional - i've tried to note where this will be, because they're not always optional now
# we'll want to be consistent with what arguments are called across functions (same name = same input)
# you'll have to put the package name in front of functions that come from a different package

# compare pairwise snv distances returned by ape::dist.dna() using pairwise.deletion = T and F
# example of using dist.dna: dist.dna(dna_aln, model = 'N', pairwise.deletion = T, as.matrix = T)

# get intra- and inter-facility pairwise snv distances (might want to change this name to explain better what we're outputting)
# dists - snv distance matrix returned by dist.dna
# should we run dist.dna in here? I guess it depends whether any other functions need that as input, then it would be in its own prep function?
# locs - locations of isolates (e.g. facility of isolation)
# pt - patient isolate was taken from (optional; will remove pairwise snv distances between the same patient)





#' Get matrix of inter- and intra-facility pariwise snv distances
#'
#' @param dists a SNV distance matrix returned by the dist.dna function from the ape package
#' @param locs a named vector of locations of isolates (e.g. facility of isolation), with the name being the sample ID
#' @param pt a named vector of patient that isolate was taken from with the name being sample ID (optional)
#'
#' @return snp_facility_paris, a matrix of inter- and intra-facility pariwise snv distances
#' @export
#'
#' @examples
get_snv_dists <- function(dists, locs, pt = NULL){
  #checks
  check_get_snv_dists_input(dists, locs, pt)

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
    snp_facility_pairs <- bind_cols(snps %>% filter(Patient1 != Patient2) %>% mutate(Pair_Type=ifelse(Loc1==Loc2,'Intra-facility pair','Inter-facility pair')))
  }
  else{
    #add labels
    snp_facility_pairs <- bind_cols(snps %>% mutate(Pair_Type=ifelse(Loc1==Loc2,'Intra-facility pair','Inter-facility pair')))
  }
  #return snp matrix
  return(snp_facility_pairs)
}

#how to plot this?
#tests both with and without pt
