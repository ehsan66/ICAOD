#  roxygen
# print ICA object
#'
#' Print method for an object of class "ICA"
#' @param x an object of class "ICA".
#' @param iter iteration.
#' @param ... argument with no further use.
#' @seealso \code{\link{mica}}, \code{\link{multica_4pl}} and \code{\link{ave}}
#' @export

print.ICA <- function(x, iter = NULL, ...){

  if (any(class(x) != c("list", "ICA")))
    stop("'x' must be of class 'ICA'")
  ## to not get confused with design points
  object <- x
  if (is.null(iter))
    totaliter <- length(object$evol) else
      totaliter <- iter

  if (totaliter > length(x$evol))
    stop("'iter' is larger than the maximum number of iterations")


  if( grepl("on_average", x$arg$type))
    type <- "optim_on_average" else
      type <- x$arg$type
  ##############################################################################
  ### printing, match with cat in iterate.ICA
  cat("\n###################################################\n## ICA generated design:")
  if (type != "locally" && type != "optim_on_average" && type != "multiple_locally")
    cat("\nICA iter:", totaliter, "\npoints:", object$evol[[totaliter]]$x,
        "\nweights: ", object$evol[[totaliter]]$w,
        "\nparam: ", object$evol[[totaliter]]$param,
        "\nbest criterion value: ", object$evol[[totaliter]]$min_cost) else
          cat("\nICA iter:", totaliter, "\nPoints:", object$evol[[totaliter]]$x,
              "\nWeights: ", object$evol[[totaliter]]$w,
              "\nBest criterion value: ", object$evol[[totaliter]]$min_cost)



  # cat("\nanswering: ")
  # cat_mat(object$evol[[totaliter]]$answering)
  # cat("\nanswering cost value: ", object$evol[[totaliter]]$answering_cost)

  cat("\nmaximum of derivative plot:", object$evol[[totaliter]]$max_deriv,
      "\nefficiency lower bound (ELB):", object$evol[[totaliter]]$ELB)

  cat("\ntotal local search:", object$alg$nlocal, "\ntotal revolution:", object$alg$nrevol)
  if (object$arg$control$only_improve)
    cat("\nTotal improve:", object$alg$nimprove)

  cat("\n###################################################")
  ##############################################################################
  return(invisible(NULL))
}

