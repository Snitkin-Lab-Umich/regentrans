#utils
#' dplyr pipe
#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

## make R CMD CHECK shut up about the dot `.``; this alone isn't working so also added names
## See: \url{https://github.com/tidyverse/magrittr/issues/29}
utils::globalVariables(c(".","subtr_size","index","isolate_name","pair_type","pairwise_dist",
                         "n","Intra-facility pair","Inter-facility pair","frac_intra",
                         "frac_inter","source_facil","n_transfers","freq","Freq","dest_facil",
                         "loc1","loc2","pairwise_dist","Patient1","Patient2","isolate1",
                         "isolate2","n_closely_related_pairs", "n_transfers_f12",
                         "n_transfers_f21", "pt_trans_metric",
                         "pt_trans_metric_f12", "pt_trans_metric_f21","isolate_id"))


#' Make matrix long-form where each row represents a cell in the matrix
#'
#' @param facil_dist symmetric or asymmetric matrix that you want to convert to long form.
#' @param col_names Column names for the output data frame (3 columns - rownames, colnames, values)
#'
#' @return long form data where each row represents a cell in the matrix (rowname, column name, value)
#' @export
#' @description Each row represents a cell in the matrix (rowname, column name, value).
#' If the matrix is symmetric, there will be one row for each symmetric pair (i.e. the number of rows in the long-form data frame will be half the number of elements in the input matrix).
#' If the matrix is asymmetric, all rows will be kept (i.e. the number of rows in the long-form data frame will be equal to the number of elements in the input matrix).
#'
#' @examples
#' make_long_form(fsp)
make_long_form <- function(facil_dist, col_names = c('loc1', 'loc2', 'fsp')){
  #check that it is a symmetric matrix
  check_long_form_input(facil_dist, col_names)

  sym <- isSymmetric(as.matrix(facil_dist))

  #change to longform
  facil_dist_long <- stats::na.omit(data.frame(as.table(as.matrix(facil_dist)))) %>% dplyr::filter(Freq != 0)
  colnames(facil_dist_long) <- col_names

  if(sym){
    ## sort facilities before summarizing (should probably make this a function)
    facil_pairs <- lapply(1:nrow(facil_dist_long), function(x)
      sort(c(as.character(facil_dist_long$loc1[x]), as.character(facil_dist_long$loc2[x])))
    )

    facil_dist_long$loc1 <- sapply(facil_pairs, function(x) x[1])
    facil_dist_long$loc2 <- sapply(facil_pairs, function(x) x[2])

    facil_dist_long <- subset_pairs(facil_dist_long)
  }

  return(facil_dist_long)
}


#' Reverse list structure
#'
#' @param ls list you want to reverse
#' @importFrom rlang :=
#'
#' @return reversed list
#' @noRd
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
