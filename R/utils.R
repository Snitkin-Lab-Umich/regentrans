#utils
#' dplyr pipe
#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

## make R CMD CHECK shut up about the dot `.``
## See: \url{https://github.com/tidyverse/magrittr/issues/29}
utils::globalVariables(c(".","subtr_size","index","isolate_name","Pair_Type","Pairwise_Dists",
                         "n","Intra-facility pair","Inter-facility pair","Frac_Intra",
                         "Frac_Inter","source_facil","n_transfers","Freq","dest_facil",
                         "Loc1","Loc2","Pairwise_Dists","Patient1","Patient2","Isolate1",
                         "Isolate2","n_closely_related_pairs", "n_transfers_f12",
                         "n_transfers_f21", "pt_trans_metric",
                         "pt_trans_metric_f12", "pt_trans_metric_f21"))

#' Reverse list structure
#'
#' @param ls list you want to reverse
#'
#' @return reversed list
#' @export
#'
#' @details Reference with example: https://stackoverflow.com/questions/15263146/revert-list-structure
reverse_list_str <- function(ls) { # @Josh O'Brien
  #checks
  check_reverse_list_str_input(ls)

  # get sub-elements in same order
  #x <- future.apply::future_lapply(ls, `[`, names(ls[[1]]))
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  #future.apply::future_apply(do.call(rbind, x), 2, as.list)
  apply(do.call(rbind, x), 2, as.list)
}
