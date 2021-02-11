# get fraction of intra-facility pairs for different snv thresholds
# dists - output of get_snv_dists - will want to change either the name of this argument or the other one so the arguments always mean the same input
# threshs - what snv thresholds to use #some vector of numbers, max number isn't > max snv distance or negative
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

#might want to change the name/call the get_snv_dists
#within it just so we don't have to run them in order,
#could write a wrapper function that calls both
