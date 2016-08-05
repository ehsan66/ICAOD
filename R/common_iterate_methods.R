#' update the object of class  'ICA'.
#'
#' Runs the algorithm for more number of iterations.
#' @param object object of class 'ICA'.
#' @param iter number of iterations.
#' @return object of class 'ICA'.
#'

iterate <- function(object, iter){
  UseMethod("iterate")
}

