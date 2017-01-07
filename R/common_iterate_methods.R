#' update the object by running the ICA algorithm for more number of iterations.
#'
#' Runs the algorithm for more number of iterations and update the results.
#' @param object object of class 'ICA'.
#' @param iter number of iterations.
#' @return object of class 'ICA'.
#' @export
#' @seealso \code{\link{iterate.ICA}}

iterate <- function(object, iter){
  UseMethod("iterate")
}

