#utils
#' dplyr pipe
#' @importFrom dplyr %>%
#' @noRd
dplyr::`%>%`

## make R CMD CHECK shut up about the dot `.``
## See: \url{https://github.com/tidyverse/magrittr/issues/29}
utils::globalVariables(c("."))

#' Reverse list structure
#'
#' @param ls list you want to reverse
#'
#' @return reversed list
#' @noRd
#'
#' @examples
#' #Reference with example: https://stackoverflow.com/questions/15263146/revert-list-structure
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
