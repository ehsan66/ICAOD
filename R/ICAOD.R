#' ICAOD: finds optimal designs for nonlinear models
#'
#' Three functions are available to find optimal designs for nonlinear models:
#' \itemize{
#'  \item{\code{\link{mica}}: }{finds locally D-optimal, minimax D-optimal and standardized maximin D-optimal designs.}
#'  \item{\code{\link{ave}}: }{finds optimum-on-the-average designs (optimal designs in average).}
#'   \item{\code{\link{multica_4pl}}: }{finds multiple-objective optimal designs for the 4-parameter logistic model or equivalently, the 4-parameter Hill model.}
#' }
#' In addition, one can use the three following functions to check the optimality of a given design using the general equivalence theorem:
#' \itemize{
#'  \item{\code{\link{equivalence}}: }{checks the optimality of a given design with respect to  the locally D-optimal, the minimax D-optimal and the standardized maximin D-optimal criteria.}
#'  \item{\code{\link{equivalence_ave}}: }{checks the optimality of a given design with respect to the optimum-on-the-average criterion.}
#'   \item{\code{\link{equivalence_multiple}}: }{checks the optimality of a given design with respect to the multiple-objective optimal criterion for the 4-parameter logistic model or equivalently, the 4-parameter Hill model.}
#' }
#'
#' @section Fisher information matrix:
#' User has to provide the Fisher information matrix for the model of interest through an argument called \code{fimfunc}.
#' This argument can be either a user-defined R function that returns the Fisher information as a matrix or a character strings that is the name of the
#' Fisher information matrix available in the ICAOD package for some popular models. See "Details" in \code{\link{mica}}.
#' @docType package
#' @name ICAOD
NULL
