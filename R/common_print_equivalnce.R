#  roxygen
#' print equivalence object
#'
#' Print method for an object of class "equivalence".
#' @param x an object of class 'equivalence'.
#' @param ... argument with no further use.
#' @seealso \code{\link{equivalence}}, \code{\link{equivalence_multiple}} and \code{\link{equivalence_on_average}}.
#' @export

print.equivalence <- function(x,  ...){

  if (any(class(x) != c("list", "equivalence")))
    stop("'x' must be of class 'equivalence'")
  ## to not get confused with design points
  object <- x
  digits1 <- 5
  ##############################################################################
  ### printing, match with cat in iterate.equivalence
  cat("\n###################################################\n## Checking equivalence theorem:")
  if (object$type == "minimax" || object$type == "standardized")
    cat("\nAll optima of the inner problem:", paste_mat(round(object$all_optima, digits1)),
        "\nCost values of all optima: ", object$all_optima_cost,
        "\nAnswerign set (chosen from all optima): ", paste_mat(round(object$answering, digits1)),
        "\nCost values of answering set: ", object$answering_cost,
        "\nProbability measure on answering set: ", object$mu)
    cat("\nCriterion value", object$crtval,
      "\nMaximum of sensitivity function: ", object$max_deriv,
      "\nD-efficiency lower bound: ", object$DLB,
      "\n###################################################")
  ##############################################################################
  return(invisible(NULL))
}

