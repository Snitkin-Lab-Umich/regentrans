# functions to include in the regional genomic transmission analysis package

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
# locs - locations of isolates (e.g. facility of isolation), named with the same names as dists
# pt - patient isolate was taken from (optional; will remove pairwise snv distances between the same patient), named with the same name as dists
get_snv_dists <- function(dists, locs, pt = NULL){
  # put checks at beginning
  check_names(rownames(dists), colnames(dists)) # see how many overlap - in this function you'll  put warnings
  # possibly make a subset function (so all are the same)
  loc_sub <- locs[rownames(dists)]
  pt_sub <- pt[rownames(dists)]
  # snps <- dists[lower.tri(dists, diag = FALSE)] <- NA
  snps <- na.omit(data.frame(as.table(as.matrix(dists))))
  snps$loc1 <- loc_sub[snps$Var1]
  snps$loc2 <- loc_sub[snps$Var2]
  snps$pt1 <- pt_sub[snps$Var1]
  snps$pt2 <- pt_sub[snps$Var2]
  snps
  bind_cols(snps %>% filter(pt1 != pt2) %>% mutate(intra=ifelse(loc1==loc2,'Intra-facility pair','Inter-facility pair')))
}

# get fraction of intra-facility pairs for different snv thresholds
# dists - output of get_snv_dists - will want to change either the name of this argument or the other one so the arguments always mean the same input
# threshs - what snv thresholds to use
get_frac_intra <- function(dists, threshs = seq(1,50,1)){
  intra_cts <- t(sapply(threshs, function(i){
    intra <- dists$intra[dists$Freq < i & !is.na(dists$intra)]
    im <- c(i,
            table(factor(intra,levels=c('Intra-facility pair','Inter-facility pair'))),
            round(mean(intra == 'Intra-facility pair'),2),
            round(mean(intra == 'Inter-facility pair'),2))
    names(im) <- c('thresh','inter','intra','frac_intra','frac_inter')
    im
  }))
  data.frame(intra_cts)
}

# get clusters
# tr - tree
# locs - locations - named vector where names are same as tree tip labels # maybe make this more general?
# pureness - how pure the cluster is (<= 0.5)
# maybe add bootstrap value
get_clusters <- function(tr, locs, pureness = 1){ # pureness shouldn't be <= 0.5?
  locs_sub <- locs[tr$tip.label]
  subtrs_sub <- ape::subtrees(tr)
  pure_subtrees <- get_largest_subtree(subtrs_sub, locs_sub, bootstrap = NULL, pureness = pureness) # NOTE: this function is in snitkitr right now, but I think we should migrate it to this package (or at least include it here as well); this _might_ be buggy, so definitely good to add unit tests for it
  pure_subtr_info <- bind_cols(f_id=locs_sub,
                               subtr_size=unlist(pure_subtrees$largest_st),
                               index=unlist(pure_subtrees$largest_st_i))
  # change singletons from 0 to 1
  pure_subtr_info <- pure_subtr_info %>% mutate(subtr_size=ifelse(subtr_size==0 & index == 1, 1, subtr_size))
  # remove duplicates (singletons aren't duplicates)
  pure_subtr_info <- pure_subtr_info[!duplicated(pure_subtr_info) | pure_subtr_info$subtr_size == 1,]
  return(pure_subtr_info=pure_subtr_info) #maybe add pureness of cluster? 
}

