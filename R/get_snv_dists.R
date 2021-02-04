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
# locs - locations of isolates (e.g. facility of isolation)
# pt - patient isolate was taken from (optional; will remove pairwise snv distances between the same patient)
get_snv_dists <- function(dists, locs, pt){
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
