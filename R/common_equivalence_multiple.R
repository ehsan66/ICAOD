#' Checking the optimality of a given design with respect to  the multi-objective criterion for the 4-parameter logisitic model.
#'
#' The equivalence theorem states
#'  that for a given vector of weights
#'  \eqn{\lambda = (\lambda_1, \lambda_2, \lambda_3)}{\lambda = (\lambda1, \lambda2, \lambda3)}, the design \eqn{\xi_\lambda}
#'  is the multi-objective optimal design if and only if for all does \eqn{x} in the dose range \eqn{\chi} (design space)
#'  \deqn{d(x, \xi_\lambda) \leq 0}{d(x, \xi_\lambda) <= 0}
#'  with equality when \eqn{x} is a dose level of design \eqn{\xi_\lambda}.
#'  See Eq. 6 of Hyun and Wong (2015) for the details.
#'
#'
#'
#' @param x a vector of design points. When design space is multi-dimensional then \code{x} should be filled dimension by dimension. See "Examples".
#' @param w a vector of design weights.
#' @param lx lower bound of the design space \eqn{\chi}.
#' @param ux upper bound of the design space \eqn{\chi}.
#' @param param initial guess for parameters \eqn{\Theta = (\theta_1, \theta_2, \theta_3, \theta_4)}{\Theta = (\theta1, \theta2, \theta3, \theta4)}.
#' @param lambda user select weights, where \eqn{\lambda_1}{\lambda1} is the weight for estimating parameters,
#' \eqn{\lambda_2}{\lambda2} is the weignt for estimating median effective dose level (ED50), and \eqn{\lambda_3}{\lambda3} is the weight for estimating minimum effective dose level (MED).
#' @param delta numeric, predetermined clinically significant effect to define the MED.
#' @param maxeval_equivalence maximum number of evaulations (\code{maxeval})  that will be passed to optimization function \code{\link[nloptr]{directL}} to find the maximum of the sensitivity function required for calculating DLB. See "Details" of \code{\link{equivalence}}.
#' @param plot_sensitivity logical; sensitivity should be plotted? see "Details" of \code{\link{equivalence}}.
#' @seealso \code{\link{equivalence}} and \code{\link{equivalence_on_average}}.
#' @details
#' When \eqn{\lambda_2 = 1}{\lambda1 = 1}, the function checks the equivalence theorem with respect to the
#'  c-optimality criterion for estimating ED50. When \eqn{\lambda_3 = 1}{\lambda3 = 1}, it checks the equivalence theorem
#'  with respect to the c-optimality criterion for estimating MED.  In both cases,
#'   due to the tolerance issue for computing the generalized inverse, the results may not be true.
#'   Therefore, this function should only be used for multiple-objective optimal design
#'   \eqn{\lambda_1 \neq 0} and  \eqn{\lambda_2 \neq 0.}{\lambda1 and \lambda2 are not equal to 0.}
#'
#'  The tolerance for finding the general inverse is set to \code{.Machine$double.xmin}.
#' @examples
#' ## verfying the design in Table 2 of Hyun and Wong (2015), first row, fisrt column.
#'Theta1 <- c(1.563, 1.790, 8.442, 0.137)
#'equivalence_multiple (x = c(log(.001), -5.21, -4.08, log(1000)),
#'                      w = c(.25, .25, .25, .25),
#'                      lx = log(.001), ux = log(1000),
#'                      param = Theta1,
#'                      lambda = c(1, 0, 0),
#'                      delta = -1)
#'
#' \dontshow{
#' \dontrun{
#'#########################################################################################
#'## examples fof using this function for c-optimal designs
#'# first row second column: c-optimal design for estimating ED50
#'equivalence_multiple (x = c(log(.001), -4.80, log(1000)),
#'                      w = c(.276, .500, .224),
#'                      lx = log(.001), ux = log(1000),
#'                      param = Theta1,
#'                      lambda = c(0, 1, 0),
#'                      delta = -1)
#'## criterion value is 1e+24 which will be returned when variance for estimating ED50 is comutationaly negative!
#'## if we change the tolerance for finding  Moore-Penrose Matrix Inverse to .Machine$double.eps
#'# when get 2.201179 for the criterion value
#'
#'equivalence_multiple (x = c(-6.910, -4.6150000, -4.6000000, 6.910),
#'                      w =   c(0.499998, 0.2230491, 0.2305728, 0.04637965),
#'                      lx = log(.001), ux = log(1000),
#'                      param = Theta1,
#'                      lambda = c(0, 0, 1),
#'                      delta = -1)
#'
#'## now let set the real value of the smallest and the largest design point! why?
#'equivalence_multiple (x = c(log(.001), -4.6150000, -4.6000000, log(1000)),
#'                      w =   c(0.499998, 0.2230491, 0.2305728, 0.04637965),
#'                      lx = log(.001), ux = log(1000),
#'                      param = Theta1,
#'                      lambda = c(0, 0, 1),
#'                      delta = -1)
#'}
#'}
#' @return
#'  an object of class \code{'equivalence'} that is a list contains:
#'  \describe{
#'  \item{\code{max_deriv}}{maximum of the sensitivity function}
#'  \item{\code{DLB}}{D-efficiency lower bound. If negative, the value of \code{maxeval_equivalence} should be increased to find the global maximum.}
#'  \item{\code{crtval}}{criterion value.}
#'  }
#' @export
equivalence_multiple <- function(x, w,
                                 lx, ux,
                                 param,
                                 lambda,
                                 delta,
                                 maxeval_equivalence = 6000,
                                 plot_sensitivity = TRUE){



  ################################
  ### fimfunc
  fimfunc <- "FIM_logistic_4par"
  lp <- up <- param
  type <- "multiple_locally"
  multiple <- list(delta = delta, lambda = lambda)
  #################################

  if (!is.function(fimfunc) && !is.character(fimfunc))
    stop("'fimfunc' can be either character or function")
  if (length(lx) != length(ux))
    stop("Length of \"lx\" is not equal to length of \"ux\"")
  if (missing(param))
    stop("\"param\" is missing")
  if (missing(lambda))
    stop("\"lambda\" is missing")



  #############################################################################
  # fimchar
  # fimfunc can be character or function. however for now it is only character.
  #if character, set fimchar and find the appropriate FIM function.
  if (is.character(fimfunc)){
    fimchar <- fimfunc
    # change 'fimfunc' to be a function corresponding to 'fimmchar'.
    fimfunc <- check_lp_fimchar(fimchar = fimchar, lp = param)
  }else
    fimchar <- "all"
  #############################################################################

  fimfunc2 <- function(x, w, param)
    fimfunc(x, w, param)
  #fimfunc(x, w, param, ...)


  temp <- create_multiple (model = fimchar, fimfunc2 = fimfunc2,  type = type, multiple = multiple)
  Psi_x <- temp$PsiMulti_x
  Psi_Mu <- temp$PsiMulti_Mu
  answering <- matrix(param, nrow = 1)
  mu <- 1

  ### finding the value of criteria


  crtval <- temp$crfunc(param = param, n_independent = 1, q = c(x, w))


  ################################################################################
  ### Psi as a function of x and x, y for plotting. Psi_x defined as minus psi to find the minimum
  ## Psi_x is mult-dimensional, x can be of two dimesnion.
  Psi_x_minus <- function(x1, mu, FIM,  x, w,  answering){
    Out <- Psi_x(x1 = x1, mu = mu, FIM = FIM,  x = x, w = w, answering = answering)
    return(-Out)
  }

  ### for plotting
  Psi_x_plot <-  Psi_x ## for PlotPsi_x


  ##########################################################################
  # find the maximum of derivative function
  OptimalityCheck <- directL(fn = Psi_x_minus,
                             lower = lx,
                             upper = ux,
                             mu = mu,
                             answering = answering,
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
                                   mu = mu,
                                   answering = answering,
                                   x = x,
                                   w = w,
                                   FIM = fimfunc2)
  check_vertices <- check_vertices$minima[, dim(check_vertices$minima)[2]]
  ## minus because for optimality check we minimize the minus derivative function
  max_deriv <- c(-OptimalityCheck$value, check_vertices)
  max_deriv <- max(max_deriv)
  ##########################################################################

  # D-efficiency lower bound
  npar <- length(param)
  DLB <- npar/(npar + max_deriv)

  if (plot_sensitivity)
    PlotPsi_x(lower = lx, upper =   ux,
              Psi_x = Psi_x_plot,
              FIM  = fimfunc2,
              mu = mu,
              x = x,
              w = w,
              plot_3d = NULL, # not applicable
              answering = answering)

  object <- list(type = "locally", max_deriv = max_deriv, DLB = DLB, crtval = crtval)
  class(object) <- c("list", "equivalence")

  return(object)
}
