######################################################################################################*
######################################################################################################*
#' ICA Control Optimization Parameters
#'
#' The function \code{ICA.control} returns a list of ICA control parameters.
#'
#'@param ncount Number of countries. Defaults to \code{40}.
#'@param nimp Number of imperialists. Defaults to 10 percent of \code{ncount}.
#'@param assim_coeff Assimilation coefficient. Defaults to \code{4}.
#'@param revol_rate Revolution rate. Defaults to \code{0.3}.
#'@param damp Damp ratio for revolution rate.  \code{revol_rate} is decreased in every iteration by a factor of \code{damp} (\code{revol_rate * damp}). Defaults to \code{0.99}.
#'@param uniting_threshold If the distance between two imperialists is less than the product of the uniting threshold by the largest distance in the search space, ICA unites the empires. Defaults to \code{0.02}.
#' @param equal_weight Should the weights of design points assumed to be equal? Defaults to \code{FALSE}. If \code{TRUE}, it reduces the dimension of the search space and produces a design that gives equal weight to all of its support points.
#' @param sym   Should the design points be symmetric around \code{sym_point}? Defaults to \code{FALSE}. When \code{TRUE}, \code{sym_point} must be given.
#' @param sym_point  If \code{sym = TRUE}, the design points will be symmetric around \code{sym_point}. See 'Details'.
#' @param stop_rule  Either  \code{'maxiter'} or \code{'equivalence'}.
#'  Denotes the type of stopping rule.  See 'Details'. Defaults to \code{'maxiter'}.
#' @param stoptol If \code{stop_rule = 'equivalence'}, algorithm stops when  ELB is larger than  \code{stoptol}. Defaults to \code{0.99}.
#' @param checkfreq The algorithm verifies the  general equivalence theorem in
#'       every \code{checkfreq} iterations.
#'       When \code{checkfreq = 0}, no verification will be done. When \code{checkfreq = Inf}, only the output design will be verified.
#'       Defaults to \code{0}.
#' @param plot_cost Plot the iterations (evolution) of algorithm? Defaults to \code{TRUE}.
#' @param plot_sens  Plot the sensitivity (derivative) function at every \code{checkfreq}. Defaults to \code{TRUE}.
#' @param plot_3d Character. Which package should be used to plot the sensitivity plot for models with two explanatory variables?
#' @param trace Print the information in every iteration? Defaults to \code{TRUE}.
#' @param rseed Random seed. Defaults to \code{NULL}.
#'@return  A list of ICA control parameters.
#' @details
#' If \code{stop_rule = 'maxiter'}, the algorithm iterates until maximum number of iterations.\cr
#'   If \code{stope_rule = 'equivalence'}, the algorithm stops when either ELB  is greater than \code{stoptol} or it reaches \code{maxiter}.
#'   In this case, you must specify the check frequency by \code{checkfreq}.
#'   Note that checking equivalence theorem is a very time consuming process,
#'    especially for Bayesian and minimax design problems.
#'   We advise using this option only for locally, multiple objective and robust optimal designs.
#'
#'
#' What to follows shows  how \code{sym_point} and \code{sym} may be useful? \cr
#'  Assume the 2PL model of the form \eqn{ P(Y=1) = \frac{1}{1+exp(-b(x - a))}}{P(Y=1) = 1/(1+exp(-b(x - a)))} and
#'  let the parameters \eqn{a} and \eqn{b}
#'   belong to
#'   \eqn{[a_L, a_U]} and \eqn{[b_L, b_U]}, respectively.
#'   It can be shown that the optimal design for this model
#'   is symmetric around \eqn{a_M = \frac{a_L + a_U}{2}}{a_M= (a_L + a_U)/2}.
#'   For this model, to find accurate symmetric designs, one can set \code{sym = TRUE} and
#'    provide the value of the \eqn{a_M} via \code{sym_point}.
#'    In this case, the output design will be symmetric around the point \code{sym_point}.
#'   The length of  \code{sym_point} must be equal to the number of model predictors, here, is equal to \code{1}.
#'
#'
#'
#'
#'@export
#'@examples
#' ICA.control(ncount = 100)
ICA.control <- function(ncount = 40, nimp = ncount/10, assim_coeff = 4, revol_rate = .3, damp = .99, uniting_threshold = .02,
                        equal_weight = FALSE, sym = FALSE, sym_point = NULL,
                        stop_rule = c("maxiter", "equivalence"), stoptol = .99, checkfreq = 0,
                        plot_cost = TRUE, plot_sens = TRUE, plot_3d = c("lattice", "rgl"),
                        trace = TRUE, rseed = NULL){

  #############################################################################*
  #functype = c("minimax", "locally", "Bayesian")
  ## ICA tuning parameters
  if (!is.numeric(ncount) || ncount <= 0)
    stop("value of 'ncount' must be > 0")
  if (!is.numeric(nimp) || nimp <= 0)
    stop("value of 'nimp' must be > 0")
  if (ncount - nimp <= nimp)
    stop( "number of colonies is less than 'nimp'. Increase 'ncount' or decrease 'nimp'")
  if (!is.numeric(assim_coeff) || assim_coeff <= 0)
    stop("value of 'assim_coeff' must be > 0")
  if (!is.numeric(revol_rate) || revol_rate <= 0)
    stop("value of 'revol_rate' must be > 0")
  if (!is.numeric(damp) || damp <= 0)
    stop("value of 'damp' must be > 0")
  #if (!is.numeric(zeta) || zeta <= 0)
  #  stop("value of 'zeta' must be > 0")
  if (!is.numeric(uniting_threshold) || uniting_threshold <= 0)
    stop("value of 'uniting_threshold' must be > 0")
  #############################################################*
  if (!is.logical(equal_weight))
    stop("'equal_weight' must be logical")
  if (!is.logical(sym))
    stop("'sym' must be logical")
  ## if TRUE, then ld and ud does not have the lower bound and upper bound for the weights.
  ## In this case, the countries are only the points and not weights
  if (sym && is.null(sym_point))
    stop("symetric point should be given by 'sym_point'")
  if (equal_weight && sym)
    stop("symmetric property does not work when only equal-weighted design is requested")
  #############################################################*
  # some functionality for stopping rules
  if(!(stop_rule[1] %in% c("maxiter", "equivalence")))
    stop("'stop_rule' can be 'maxiter' or 'equivalence'")
  # if (functype[1] == "locally")
  #   stop_rule <- "equivalence"
  if (!is.numeric(stoptol) || stoptol <= 0 || stoptol > 1)
    stop("value of 'stoptol' must be > 0 and <= 1")
  if (!is.numeric(checkfreq) || checkfreq < 0)
    stop("value of 'checkfreq' must be >= 0 and a numeric")
  # if (functype[1] == "Bayesian")
  #   checkfreq <- Inf
  if (plot_sens && checkfreq == 0)
    warnings("the sensitivity plot will not be plotted when 'checkfreq' in 'ICA.control' is zero.")
  if (!all(plot_3d %in%  c("lattice", "rgl")))
    stop("'plot_3d' must be 'lattice' or 'rgl'")
  ############################################################*

  if (!is.logical(plot_cost))
    stop("'plot_cost' must be logical")
  # if (functype[1] == "locally")
  #   plot_cost <- FALSE
  if (!is.logical(plot_sens))
    stop("'plot_sens' must be logical")
  if (!is.logical(trace))
    stop("'trace' must be logical")
  # if (functype[1] == "locally")
  #   trace <- FALSE
  return(list(ncount = ncount, nimp = nimp, assim_coeff = assim_coeff, revol_rate = revol_rate, damp = damp,  uniting_threshold = uniting_threshold,
              equal_weight = equal_weight, sym = sym, sym_point = sym_point,
              stop_rule = stop_rule[1], stoptol = stoptol, checkfreq = checkfreq,
              plot_cost = plot_cost, plot_sens = plot_sens, trace = trace, rseed = rseed))
}
######################################################################################################*
######################################################################################################*

######################################################################################################*
######################################################################################################*
#'@title Updating an Object of Class \code{'bayes'} or \code{'minimax'}
#' @description  Runs ICA for more number of iterations.
#' @param object An object of class \code{'minimax'} or  \code{'bayes'}.
#' @param iter Number of iterations.
#' @return An (updated) object of class \code{'bayes'} or \code{'minimax'}.
#' @export
#' @seealso \code{\link{iterate.minimax}} and \code{\link{iterate.bayes}}.
iterate <- function(object, iter){
  UseMethod("iterate")
}

