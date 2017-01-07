#' Checking the optimality of a  design with respect to the optim-on-the-average criterion
#'
#' Let \eqn{\Theta}  be the set of plausible parameter values and weighted by a probability measure
#' \eqn{\pi}, the measure having support in the parameter space \eqn{\Theta}, and \eqn{\Psi(\xi, \theta) =|M(\xi, \theta)|}.
#' A design \eqn{\xi^*}{\xi*} is optimum-on-the-average with respect to prior \eqn{\pi}  if
#'  the following inequality holds for all \eqn{\boldsymbol{x} \in \chi}{x belong to \chi}
#'  \deqn{c(\boldsymbol{x}, \pi, \xi^*) = \int_{\pi} tr M^{-1}(\xi^*, \theta)I(\boldsymbol{x}, \theta)\pi(\theta) d(\theta)-p \leq 0,}{
#'          c(x, \pi, \xi*)  = integration over \pi with integrand tr M^-1(\xi*, \theta)I(x, \theta)\pi(\theta) d(\theta)-p <= 0}
#'           with equality at all support points of \eqn{\xi^*}{\xi*}.
#'            Here, \eqn{p} is the number of model parameters.
#'
#'
#' @param fimfunc Fisher information matrix. Can be the name of the Fisher information matrix from FIM family functions available in this package as a
#'  character string or a function that returns the information matrix. See "Details" of \code{\link{mica}}.
#' @param x a vector of design points. If the  model has \eqn{n} explanatory variables, let \eqn{x_{ij}}
#'  be the \eqn{j}th component of the $\eqn{i}th design point.
#' The argument \code{x} is \eqn{x = (x_{11}, x_{21},..., x_{k1},..., x_{1n}, x_{2n},... x_{kn})}.
#' See "Examples" on how to set this argument when the design space does not have one dimension, e.g. is of two-dimension.
#' @param w a vector of design weights
#' @param lx lower bound of the design space \eqn{\chi}
#' @param ux upper bound of the design space \eqn{\chi}
#' @param prior a vector of the probability measure \eqn{\pi}.
#' @param param a matrix for set of parameters, i.e. support of \eqn{\pi}. Every row is is a vector of values of a parameter.
#' The number of its rows must be equal to the length of \code{prior}.
#' @param maxeval_equivalence maximum number of evaluations (\code{maxeval})  that will be passed to optimization
#' function \code{\link[nloptr]{directL}} to find the maximum of the sensitivity function required for calculating ELB.
#'  See "Details" of \code{\link{equivalence}}.
#' @param plot_sensitivity  \code{logical}, if \code{TRUE}, the sensitivity function will be plotted.
#' @param ... further argument to be passed to \code{fimfunc}
#' @examples
#'equivalence_ave(fimfunc ="FIM_logistic",lx = -5, ux = 5, x = c(0.2603688, 1, 1.739631),
#'                       w = c(0.2750147, 0.4499705, 0.2750148),  prior = c(.25, .25, .25, .25),
#'                       param =  matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2))
#' @return
#'  an object of class \code{'equivalence'} that is a list contains:
#'  \describe{
#'  \item{\code{max_deriv}}{maximum of the sensitivity function}
#'  \item{\code{ELB}}{Efficiency lower bound. If it is negative,
#'   then the value of \code{maxeval_equivalence} should be increased to find the global maximum.}
#'  \item{\code{crtval}}{criterion value}
#'  }
#' @seealso \code{\link{equivalence}} and \code{\link{equivalence_multiple}}.
#' @export
equivalence_ave <- function(fimfunc,
                                   x, w,
                                   lx, ux,
                                   prior,
                                   param,
                                   maxeval_equivalence = 6000,
                                   plot_sensitivity = TRUE,
                                   ...){

  if (!is.function(fimfunc) && !is.character(fimfunc))
    stop("'fimfunc' can be either character or function")
  if (length(lx) != length(ux))
    stop("Length of \"lx\" is not equal to length of \"ux\"")
  if (missing(param))
    stop("\"param\" is missing")
  if (missing(prior))
    stop("\"prior\" is missing")
  if (length(prior) != dim(param)[1])
    stop("length of \"prior\" is not equal to the number of rows of \"param\"")

  #############################################################################
  # fimchar
  # fimfunc can be character or function. however for now it is only character.
  #if character, set fimchar and find the appropriate FIM function.
  if (is.character(fimfunc)){
    fimchar <- fimfunc
    # change 'fimfunc' to be a function corresponding to 'fimmchar'.
    fimfunc <- check_lp_fimchar(fimchar = fimchar, lp = param[1, ])
  }else
    fimchar <- "all"
  #############################################################################


  fimfunc2 <- function(x, w, param)
    fimfunc(x, w, param, ...)
  ###############################################################################
  ## needed for reporting the value of the criterion
    crfunc_on_average_D <- function(param, q, n_independent) {
      lq <- length(q) # q is the design points and weights
      pieces <- lq / (n_independent + 1)
      x_ind <- 1:(n_independent * pieces)
      w_ind <- (x_ind[length(x_ind)] + 1):lq
      x <- q[x_ind]
      w <- q[w_ind]
      on_average_crfunc <- sum(
        sapply(1:length(prior), FUN= function(j)
          prior[j] * -det2(fimfunc2(x = x, w = w, param = param[j, ]), logarithm = TRUE))
      )
      return(on_average_crfunc)
    }
  ###############################################################################


  ################################################################################
  ### Psi as a function of x and x, y for plotting. Psi_x defined as minus psi to find the minimum
  ## Psi_x is mult-dimensional, x can be of two dimesnion.
  Psi_x_minus <- function(x1, mu, FIM,  x, w,  answering){
    Out <- Psi_x(x1 = x1, mu = mu, FIM = FIM,  x = x, w = w, answering = answering)
    return(-Out)
  }

  ### for plotting
  if(length(lx) == 1)
    Psi_x_plot <-  Psi_x ## for PlotPsi_x
  if(length(lx) == 2)
    Psi_x_plot <- Psi_xy

  ##########################################################################
  # find the maximum of derivative function
  OptimalityCheck <- directL(fn = Psi_x_minus,
                             lower = lx,
                             upper = ux,
                             mu = prior,
                             answering = param,
                             x = x,
                             w = w,
                             FIM = fimfunc2,
                             nl.info = FALSE,
                             control=list(xtol_rel=.Machine$double.eps,
                                          maxeval = maxeval_equivalence))
  ##sometimes the optimization can not detect maximum  in the bound, so here we add the cost values on the bound
  vertices_outer <- make_vertices(lower = lx, upper = ux)
  check_vertices <- find_on_points(fn = Psi_x,
                                   points = vertices_outer,
                                   mu = prior,
                                   answering = param,
                                   x = x,
                                   w = w,
                                   FIM = fimfunc2)
  check_vertices <- check_vertices$minima[, dim(check_vertices$minima)[2]]
  ## minus because for optimality check we minimize the minus derivative function
  max_deriv <- c(-OptimalityCheck$value, check_vertices)
  max_deriv <- max(max_deriv)
  ##########################################################################

  # D-efficiency lower bound
  npar <- dim(param)[2]
  ELB <- npar/(npar + max_deriv)

  if (plot_sensitivity)
    PlotPsi_x(lower = lx, upper =   ux,
              Psi_x = Psi_x_plot,
              FIM  = fimfunc2,
              mu = prior,
              x = x,
              w = w,
              plot_3d = "lattice", # not applicable
              answering = param)


  crtval <- crfunc_on_average_D(param = param, q = c(x, w), n_independent = 1)
  object <- list(type = "locally", max_deriv = max_deriv, ELB = ELB, crtval = crtval)
  class(object) <- c("list", "equivalence")

  return(object)
}
