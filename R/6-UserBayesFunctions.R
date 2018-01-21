######################################################################################################*
######################################################################################################*
#' @title  Bayesian D-Optimal Designs
#'
#' @inheritParams minimax
#' @param prior An object of class \code{cprior}. User can also use one of the functions
#'  \code{\link{uniform}}, \code{\link{normal}},
#' \code{\link{skewnormal}} or \code{\link{student}}  to create the  prior. See 'Details' of \code{\link{bayes}}.
#' @param crt.bayes.control Control parameters to approximate the integral in  Bayesian criterion at a given design over the parameter space.
#'  For details, see \code{\link{crt.bayes.control}}.
#' @param sens.bayes.control Control parameters to verify the general equivalence theorem. For details, see \code{\link{sens.bayes.control}}.
#' @param npar Number of model parameters.  Used when \code{fimfunc} is given instead of \code{formula} to specify the number of model parameters.
#'   If not specified truly, the sensitivity (derivative) plot may be shifted below the y-axis. When \code{NULL}, it will be set to \code{length(parvars)} or
#'   \code{prior$npar} when \code{missing(formula)}.
#' @export
#' @description
#'  Finds (pseudo) Bayesian D-optimal designs for nonlinear models.
#'  It should be used when the user assumes a (truncated) prior distribution for the unknown model parameters.
#'  If the prior is discrete, please use \code{\link{robust}}.
#'
#'
#' @details
#'  Let \eqn{\Xi} be the space of all  approximate designs with
#'  \eqn{k} design points (support points) at \eqn{x_1, x_2, ...,  x_k}{x1, x2, ...,  xk} from  design space \eqn{\chi} with
#'  corresponding weights  \eqn{w_1, . . . ,w_k}{w1, . . . ,wk}.
#'  Let \eqn{M(\xi, \theta)} be the Fisher information
#'   matrix (FIM) of a \eqn{k-}point design \eqn{\xi}
#'   and  \eqn{\pi(\theta)} is a user-given  prior distribution for the vector of unknown parameters \eqn{\theta}.
#'  A  Bayesian D-optimal design \eqn{\xi^*}{\xi*} minimizes over \eqn{\Xi}
#'   \deqn{\int_{\theta \in \Theta} -\log|M(\xi, \theta)| \pi(\theta) d\theta.}{
#'    integration over \Theta -log|M(\xi, \theta)|\pi(\theta) d\theta.}
#'
#' An object of class \code{cprior}  is a  list with the following components:
#' \itemize{
#'  \item{fn: }{Prior distribution as an R \code{function} with argument \code{param} that is the vector of the unknown parameters. See below.}
#'  \item{npar: }{Number of unknown parameters and is equal to the length of \code{param}}.
#'  \item{lower: }{Argument \code{lower}. It has the same length as \code{param}}.
#'  \item{upper: }{Argument \code{lower}. It has the same length as \code{param}}.
#' }
#' A \code{cprior} object  will be passed to the argument \code{prior} of the function \code{\link{bayes}}.
#'  The argument \code{param} in \code{fn} has the same order as the argument \code{parvars} when the model is specified by a \code{formula}.
#' Otherwise, it is same as the argument \code{param} in the function \code{fimfunc}.\cr
#' The user can apply the available prior (object \code{cprior}) creators that are \code{\link{uniform}}, \code{\link{normal}},
#' \code{\link{skewnormal}} and \code{\link{student}} to create a \code{cprior} object.
#'
#' Use \code{\link{plot}} function to verify the general equivalence theorem for the output design.
#'
#' \strong{To increase the speed of the algorithm, change the tuning parameters \code{tol} and \code{maxEval} via
#' the argument  \code{crt.bayes.control}.}
#' In this case, the user should find a trade-off between accuracy and speed for his/her example.
#'
#' If some of the parameters are fixed in a model, they should be set
#' to their values via the argument \code{paravars}. In this case,
#' you must provide the number of parameters via argument \code{npar} for verifying the general equivalence theorem.
#'  See 'Examples', Section 'Weibull',  'Richards' and 'Exponential' model.
#' @return
#'  an object of class \code{bayes} that is a list including three sub-lists:
#' \describe{
#'   \item{\code{arg}}{A list of design and algorithm parameters.}
#'   \item{\code{evol}}{A list of length equal to the number of iterations that stores the information about the best design (design with least criterion value) of each iteration as follows:
#'    \code{evol[[iter]]} contains:
#'     \tabular{lll}{
#'       \code{iter}                   \tab      \tab Iteration number.\cr
#'       \code{x}                      \tab      \tab Design points. \cr
#'       \code{w}                      \tab      \tab Design weights. \cr
#'       \code{min_cost}               \tab      \tab Cost (Bayesian criterion value) of the best imperialist.  \cr
#'       \code{mean_cost}              \tab      \tab Mean of costs of all imperialists. \cr
#'       \code{sens}                   \tab      \tab An object of class 'sensbayes'. See below.\cr
#'     }
#'   }
#'
#'   \item{\code{empires}}{A list of all empires of the last iteration.}
#'   \item{\code{alg}}{A list with following information:
#'     \tabular{lll}{
#'       \code{nfeval}           \tab      \tab Number of function evaluations. See below. \cr
#'       \code{nlocal}           \tab      \tab Number of successful local search. \cr
#'       \code{nrevol}           \tab      \tab Number of successful revolutions. \cr
#'       \code{nimprove}         \tab      \tab Number of successful movements toward the imperialists in assimilation step. \cr
#'       \code{convergence}      \tab      \tab Stopped by \code{'maxiter'} or \code{'equivalence'}?\cr
#'     }
#'   }
#' }
#' \code{sens}  contains information about design verification by the general equivalence theorem.
#'  See \code{sensbayes} for more Details. It is only available every \code{ICA.control$checkfreq} iterations
#'  and the last iteration if   \code{ICA.control$checkfreq >= 0}. Otherwise, \code{NULL}.
#'
#'  \code{nfeval} does not count the function evaluations from checking the general equivalence theorem.
#'
#' @example inst/examples/bayes_examples.R
#' @seealso \code{\link{sensbayes}}
#' @importFrom cubature hcubature
#' @importFrom stats gaussian
#' @importFrom stats binomial
bayes <- function(formula,
                  predvars,
                  parvars,
                  family = gaussian(),
                  prior,
                  lx,
                  ux,
                  iter,
                  k,
                  fimfunc = NULL,
                  ICA.control =  list(),
                  crt.bayes.control = list(),
                  sens.bayes.control = list(),
                  initial = NULL,
                  npar = NULL,
                  plot_3d = c("lattice", "rgl")) {

  if (is.null(npar)){
    if (!missing(formula))
      npar <- length(parvars)
    else
      npar <- prior$npar
  }
  if (!is.numeric(npar))
    stop("'npar' is the number of parameters and must be numeric")
  output <-  bayes_inner(fimfunc = fimfunc,
                         formula = formula,
                         predvars = predvars,
                         parvars = parvars,
                         family = family,
                         lx = lx,
                         ux = ux,
                         type = "D",
                         method = "cubature",
                         iter = iter,
                         k = k,
                         npar = npar,
                         prior = prior,
                         compound = list(prob = NULL, alpha = NULL),
                         multiple.control = list(),
                         ICA.control =  ICA.control,
                         crt.bayes.control = crt.bayes.control,
                         sens.bayes.control = sens.bayes.control,
                         initial = initial,
                         plot_3d = plot_3d[1],
                         const = list(ui = NULL, ci = NULL, coef = NULL))
  return(output)
}
######################################################################################################*
######################################################################################################*
#' @title Verifying Optimality of Bayesian D-optimal Designs
#' @inheritParams bayes
#' @inheritParams sensminimax
#' @description
#'  This function plot the sensitivity (derivative) function given an approximate (continuous) design and calculate the efficiency lower bound (ELB) for Bayesian D-optimal designs.
#' Let \eqn{\boldsymbol{x}}{x} belongs to \eqn{\chi} that denotes the design space.
#' Based on the general equivalence theorem, generally, a design \eqn{\xi^*}{\xi*} is optimal if and only if the value of its sensitivity (derivative) function
#' be non-positive for all \eqn{\boldsymbol{x}}{x} in \eqn{\chi} and it only reaches zero
#' when \eqn{\boldsymbol{x}}{x} belong to the support of \eqn{\xi^*}{\xi*} (be equal to one of the design point).
#' Therefore, the user can look at the sensitivity plot and the ELB and decide whether the
#' design is optimal or close enough to the true optimal design (ELB tells us that without knowing the latter).
#'@details
#' Let \eqn{\Xi} be the space of all  approximate designs with
#'  \eqn{k} design points (support points) at \eqn{x_1, x_2, ...,  x_k}{x1, x2, ...,  xk} from  design space \eqn{\chi} with
#'  corresponding weights  \eqn{w_1, . . . ,w_k}{w1, . . . ,wk}.
#'  Let \eqn{M(\xi, \theta)} be the Fisher information
#'   matrix (FIM) of a \eqn{k-}point design \eqn{\xi}
#'   and  \eqn{\pi(\theta)} is a user-given  prior distribution for the vector of unknown parameters \eqn{\theta}.
#' A design \eqn{\xi^*}{\xi*} is Bayesian D-optimal among all designs on \eqn{\chi} if and only if  the following inequality holds for all \eqn{\boldsymbol{x} \in \chi}{x belong to \chi}
#'  \deqn{c(\boldsymbol{x}, \xi^*) =  \int_{\theta \in Theta}tr M^{-1}(\xi^*, \theta)I(\boldsymbol{x}, \theta)-p \pi(\theta) d\theta\leq 0,}{
#'  c(x, \xi*) =  integration over \Theta tr M^-1(\xi*, \theta)I(x, \theta)-p <= 0,}
#'  with equality at all support points of \eqn{\xi^*}{\xi*}.
#'  Here, \eqn{p} is the number of model parameters.
#'  \eqn{c(\boldsymbol{x},\xi^*)}{c(x, \xi*)} is
#'   called \strong{sensitivity} or \strong{derivative} function.
#'
#'  \strong{Sometimes, the CPU time can be considerably reduced
#' by choosing less conservative values for the tuning parameters \code{tol} and \code{maxEval} in
#' the function \code{\link{sens.bayes.control}}.}
#' The user should find a trade-off between accuracy and speed for his/her problem.
#'  See 'Examples'.
#' @note
#' Having accurate plots for the sensitivity (derivative) function
#'  and calculating ELB to a high precision is the primary goal here,
#'   although the process may take too long (even hours) due to
#' requesting very accurate integral approximations.
#'@export
#'@example inst/examples/sensbayes_examples.R
sensbayes <- function(formula,
                      predvars, parvars,
                      family = gaussian(),
                      x, w,
                      lx, ux,
                      fimfunc = NULL,
                      prior = list(),
                      sens.bayes.control = list(),
                      crt.bayes.control = list(),
                      plot_3d = c("lattice", "rgl"),
                      plot_sens = TRUE,
                      npar = NULL,
                      calculate_criterion = TRUE,
                      silent = FALSE){


  if (is.null(npar)){
    if (!missing(formula))
      npar <- length(parvars)
    else
      npar <- prior$npar
  }
  if (!is.numeric(npar))
    stop("'npar' is the number of parameters and must be numeric")
  output <- sensbayes_inner (formula = formula,
                             predvars = predvars, parvars = parvars,
                             family =  family,
                             x = x, w = w,
                             lx = lx, ux = ux,
                             fimfunc = fimfunc,
                             prior = prior,
                             sens.bayes.control = sens.bayes.control,
                             crt.bayes.control = crt.bayes.control,
                             type = "D",
                             plot_3d = plot_3d[1],
                             plot_sens =  plot_sens,
                             const = list(ui = NULL, ci = NULL, coef = NULL),
                             compound = list(prob = NULL, alpha = NULL),
                             varlist = list(),
                             calledfrom = "sensfuncs",
                             npar = npar,
                             calculate_criterion = calculate_criterion,
                             silent = silent)
  return(output)
}
######################################################################################################*
######################################################################################################*
#' @title   Bayesian Compound DP-Optimal Designs
#'
#' @description
#'  Finds compound Bayesian DP-optimal designs that meets the dual goal of the parameter estimation and
#'   increasing the probability of a particular outcome in a binary response  model.
#'A compound Bayesian DP-optimal design maximizes  the product of the efficiencies of a design \eqn{\xi} with respect to D- and average P-optimality, weighted by a pre-defined mixing constant
#' \eqn{0 \leq \alpha \leq 1}{0 <= \alpha <= 1}.
#'
#' @inheritParams bayes
#' @param alpha A value between 0 and 1.
#' Compound or combined DP-criterion  is the product of the efficiencies of a design  with respect to D- and average P- optimality, weighted by \code{alpha}.
#' @param prob Either \code{formula} or a \code{function}. when function, its argument is \code{x} and \code{param} same as the arguments in \code{fimfunc}.
#' \code{prob} as a function takes the design points and vector of parameters and returns the probability of success at each design points.
#' See 'Examples'.
#' @details
#' Let \eqn{\Xi} be the space of all  approximate designs with
#'  \eqn{k} design points (support points) at \eqn{x_1, x_2, ...,  x_k}
#'   from  design space \eqn{\chi} with
#'  corresponding weights  \eqn{w_1,... ,w_k}.
#'  Let \eqn{M(\xi, \theta)} be the Fisher information
#'   matrix (FIM) of a \eqn{k-}point design \eqn{\xi},
#'    \eqn{\pi(\theta)} is a user-given  prior distribution for the vector of unknown parameters \eqn{\theta} and
#'    \eqn{p(x_i, \theta)} is the ith probability of success
#' given by \eqn{x_i} in a binary response model.
#'   A Bayesian compound DP-optimal criterion maximizes over \eqn{\Xi}
#' \deqn{\int_{\theta \in \Theta} \frac{\alpha}{q}\log|M(\xi, \theta)| + (1- \alpha)
#'\log \left( \sum_{i=1}^k w_ip(x_i, \theta) \right) \pi(\theta) d\theta.}{
#' integration over \Theta \alpha/q log|M(\xi, \theta)| + (1- \alpha)
#'log ( \sum w_i p(x_i, \theta)) \pi(\theta) d\theta.
#'}
#'
#' Use \code{\link{plot}} function to verify the general equivalence theorem for the output design.
#'
#' To increase the speed of the algorithm, change the tuning parameters \code{tol} and \code{maxEval} via the
#' argument \code{crt.bayes.control}.
#'
#' @references  McGree, J. M., Eccleston, J. A., and Duffull, S. B. (2008). Compound optimal design criteria for nonlinear models. Journal of Biopharmaceutical Statistics, 18(4), 646-661.
#' @export
#' @inherit bayes return
#' @example inst/examples/bayescomp_examples.R
#' @seealso \code{\link{sensbayescomp}}
bayescomp <- function(formula,
                      predvars,
                      parvars,
                      family = binomial(),
                      prior,
                      alpha,
                      prob,
                      lx,
                      ux,
                      iter,
                      k,
                      fimfunc = NULL,
                      ICA.control =  list(),
                      crt.bayes.control = list(),
                      sens.bayes.control = list(),
                      initial = NULL,
                      npar = NULL,
                      plot_3d = c("lattice", "rgl")) {

  if (is.null(npar)){
    if (!missing(formula))
      npar <- length(parvars)
    else
      npar <- prior$npar
  }
  if (!is.numeric(npar))
    stop("'npar' is the number of parameters and must be numeric")
  if (is.formula(prob)){
    prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
  }else{
    if (!is.function(prob))
      stop("'prob' must be either a function or a formula")
    if (!formalArgs(prob) %in% c("x", "param"))
      stop("arguments of 'prob' must be 'x' and 'param'")
  }

  output <-  bayes_inner(fimfunc = fimfunc,
                         formula = formula,
                         predvars = predvars,
                         parvars = parvars,
                         family = family,
                         lx = lx,
                         ux = ux,
                         type = "DPA",
                         method = "cubature",
                         iter = iter,
                         k = k,
                         npar = npar,
                         prior = prior,
                         compound = list(prob = prob, alpha = alpha),
                         multiple.control = list(),
                         ICA.control =  ICA.control,
                         crt.bayes.control = crt.bayes.control,
                         sens.bayes.control = sens.bayes.control,
                         initial = initial,
                         const = list(ui = NULL, ci = NULL, coef = NULL),
                         plot_3d = plot_3d[1])
  return(output)
}

######################################################################################################*
######################################################################################################*
######################################################################################################*
######################################################################################################*
#'@title Verifying Optimality of Bayesian Compound DP-optimal Designs
#'@description
#'  This function plot the sensitivity (derivative) function given an approximate (continuous) design and calculate the efficiency lower bound (ELB) for Bayesian D-optimal designs.
#' Let \eqn{\boldsymbol{x}}{x} belongs to \eqn{\chi} that denotes the design space.
#' Based on the general equivalence theorem, generally, a design \eqn{\xi^*}{\xi*} is optimal if and only if the value of its sensitivity (derivative) function
#' be non-positive for all \eqn{\boldsymbol{x}}{x} in \eqn{\chi} and it only reaches zero
#' when \eqn{\boldsymbol{x}}{x} belong to the support of \eqn{\xi^*}{\xi*} (be equal to one of the design point).
#' Therefore, the user can look at the sensitivity plot and the ELB and decide whether the
#' design is optimal or close enough to the true optimal design (ELB tells us that without knowing the latter).
#'
#'@inheritParams sensbayes
#'@inheritParams bayescomp
#'@inherit sensbayes return
#'@export
#'@details
#'  \strong{Sometimes, the CPU time can be considerably reduced
#' by choosing less conservative values for the tuning parameters \code{tol} and \code{maxEval} in
#' the function \code{\link{sens.bayes.control}}.}
#' The user should find a trade-off between accuracy and speed for his/her problem.
#' @note
#' Having accurate plots for the sensitivity (derivative) function
#'  and calculating ELB to a high precision is the primary goal here,
#'   although the process may take too long (even hours) due to
#' requesting very accurate integral approximations.
#'
#'
#'
#' @seealso \code{\link{bayescomp}}
#'@example inst/examples/sensbayescomp_examples.R
sensbayescomp <- function(formula,
                          predvars, parvars,
                          family = gaussian(),
                          x, w,
                          lx, ux,
                          fimfunc = NULL,
                          prior = list(),
                          prob, alpha,
                          sens.bayes.control = list(),
                          crt.bayes.control = list(),
                          plot_3d = c("lattice", "rgl"),
                          plot_sens = TRUE,
                          npar = NULL,
                          calculate_criterion = TRUE,
                          silent = FALSE){


  if (is.null(npar)){
    if (!missing(formula))
      npar <- length(parvars)
    else
      npar <- prior$npar
  }
  if (!is.numeric(npar))
    stop("'npar' is the number of parameters and must be numeric")
  if (is.formula(prob)){
    prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
  }else{
    if (!is.function(prob))
      stop("'prob' must be either a function or a formula")
    if (!args(prob) %in% c("x", "param"))
      stop("arguments of 'prob' must be 'x' and 'param'")
  }
  output <- sensbayes_inner (formula = formula,
                             predvars = predvars, parvars = parvars,
                             family =  family,
                             x = x, w = w,
                             lx = lx, ux = ux,
                             fimfunc = fimfunc,
                             prior = prior,
                             sens.bayes.control = sens.bayes.control,
                             crt.bayes.control = crt.bayes.control,
                             type = "DPA",
                             plot_3d = plot_3d[1],
                             plot_sens =  plot_sens,
                             const = list(ui = NULL, ci = NULL, coef = NULL),
                             compound = list(prob = prob, alpha = alpha),
                             varlist = list(),
                             calledfrom = "sensfuncs",
                             npar = npar,
                             calculate_criterion = calculate_criterion,
                             silent  = silent)
  return(output)
}
######################################################################################################*
######################################################################################################*
# beff <- function(fim, xopt, wopt, x, w, prior, control = list()){
#   ### relative efficieny of x with respect to xopt
#   if (is.null(control$tol))
#     control$tol = 1e-5
#   if (is.null(control$maxEval))
#     control$maxEval = 50000
#   truncated_standard <- cubature::hcubature(f = function(param) prior$fn(t(param)),
#                                             lowerLimit = prior$lower,
#                                             upperLimit = prior$upper,
#                                             vectorInterface = TRUE)$integral
#
#
#   cr_integrand <- function(param, x, w){
#     bcrfunc1 <- apply(param, 2,
#                       FUN = function(col_par)-det2(fim(x = x, w = w, param = col_par), logarithm = TRUE)) * prior$fn(t(param))
#     dim(bcrfunc1) <- c(1, length(bcrfunc1))
#     return(bcrfunc1)
#   }
#
#
#   crfunc_bayesian_D  <- function(x, w, maxEval, tol) {
#     out <- cubature::hcubature(f = cr_integrand, lowerLimit = prior$lower, upperLimit = prior$upper,
#                                vectorInterface = TRUE,
#                                x = x, w = w, tol = control$tol, maxEval = control$maxEval)
#
#     val <- out$integral
#     return(list(val = val, fneval = out$functionEvaluations))
#   }
#
#   releff<- crfunc_bayesian_D(x = xopt, w = wopt, maxEval = control$maxEval, tol = control$tol)$val/
#     crfunc_bayesian_D(x = x, w = w, maxEval = control$maxEval, tol = control$tol)$val
#   return(releff)
# }
######################################################################################################*
######################################################################################################*
#  roxygen
#' Plotting \code{bayes} Objects
#'
#' @param x An object of class \code{bayes}.
#' @param iter Iteration number. if \code{NULL}, will be set to last iteration.
#' @param sensitivity Logical. If \code{TRUE}, the general equivalence theorem is used to check the optimality if the best design in iteration number \code{iter} and the sensitivity plot will be plotted.
#' @param calculate_criterion Re-calculate the criterion value? It only assumes a continuous parameter space for the minimax and standardized maximin designs.  Defaults to \code{FALSE}. See 'Details'.
#' @param sens.bayes.control Control parameters to verify general equivalence theorem. For details, see \code{\link{sens.bayes.control}}.
#' @param crt.bayes.control Control parameters to approximate the integration in the Bayesian criterion at a given design.
#'  For details, see \code{\link{crt.bayes.control}}.
#' @param silent Do not print anything? Defaults to \code{FALSE}.
#' @param plot_3d Which package should be used to plot the sensitivity function for two-dimensional design space. Defaults to \code{plot_3d = "lattice"}.
#' Only applicable when \code{sensitivity = TRUE}.
#' @param evolution Plot Evolution? Defaults to \code{FALSE}.
#' @param ... Argument with no further use.
#' @seealso \code{\link{bayes}}, \code{\link{bayescomp}}
#' @description
#'  This function plots the evolution of the algorithm till iteration number \code{iter} iteration and re-checks the general equivalence theorem by plotting the sensitivity function and calculating the ELB.
#' @details
#'  The criterion value can also be re-calculated for the output designs using new set of tuning parameters in the function \code{\link{crt.bayes.control}}.
#'  This is useful for  Bayesian optimal designs to assess the robustness of the
#'  criterion value with respect to different values of the tuning parameters.
#'  To put it simple, for these designs, the user can re-calculate the
#'  criterion value (approximate the integration given an output design in a Bayesian problem) with different values for  \code{maxEval} and \code{tol} in \code{\link{crt.bayes.control}}
#'  to be sure that the function \code{hcubature} approximates the integrals to an acceptable accuracy using the default values
#'  (or a new user-given values, in case it has been reset) of \code{maxEval} and \code{tol}.
#' @export
plot.bayes <- function(x, iter = NULL,
                       sensitivity = TRUE,
                       calculate_criterion = FALSE,
                       sens.bayes.control = list(),
                       crt.bayes.control = list(),
                       silent = FALSE,
                       plot_3d = c("lattice", "rgl"),
                       evolution = FALSE,
                       ...){
  if (!evolution & !sensitivity){
    warning("Both 'sensitivity' and 'evolution' set to be FALSE.\nNo action is done in plot function!")
    return(invisible(NULL))
  }
  if (any(class(x) != c("list", "bayes")))
    stop("'x' must be of class 'bayes'")
  ## to not be confused with design points
  obj <- x
  arg <- obj$arg
  if(is.null(iter))
    totaliter <- length(obj$evol) else
      totaliter <- iter
  if (totaliter > length(x$evol))
    stop("'iter' is larger than the maximum number of iterations")


  if (calculate_criterion || sensitivity){
    if (is.null(sens.bayes.control)){
      sens.bayes.control <- arg$sens.bayes.control}
    else {
      sens.bayes.control <- do.call("sens.bayes.control", sens.bayes.control)
    }
    if (is.null(crt.bayes.control)){
      crt.bayes.control <- arg$crt.bayes.control
      Psi_x_bayes  <- arg$Psi_funcs$Psi_x_bayes
      Psi_xy_bayes  <- arg$Psi_funcs$Psi_xy_bayes
    } else {

      crt.bayes.control <- do.call("crt.bayes.control", crt.bayes.control)
      temp_psi <- create_Psi_bayes(type = arg$type, prior = arg$prior, FIM = arg$FIM, lp = arg$prior$lower,
                                   up = arg$prior$upper, npar = arg$npar,
                                   truncated_standard = arg$truncated_standard,
                                   const = arg$const, sens.bayes.control = sens.bayes.control,
                                   compound = arg$compound)
      Psi_x_bayes  <- temp_psi$Psi_x_bayes
      Psi_xy_bayes  <- temp_psi$Psi_xy_bayes
    }

    vertices_outer <- make_vertices(lower = arg$lx, upper = arg$ux)
    sens_varlist <-list(npred = length(arg$lx),
                        # plot_3d = "lattice",
                        npar = arg$npar,
                        fimfunc_sens = arg$FIM_sens,
                        Psi_x_bayes  = Psi_x_bayes,
                        Psi_xy_bayes  = Psi_xy_bayes,
                        crfunc = arg$crfunc,
                        vertices_outer = vertices_outer)

    sens_res <- sensbayes_inner(x = obj$evol[[totaliter]]$x, w = obj$evol[[totaliter]]$w,
                                lx = arg$lx, ux = arg$ux,
                                fimfunc = arg$FIM,
                                prior = arg$prior,
                                sens.bayes.control = sens.bayes.control,
                                crt.bayes.control = crt.bayes.control,
                                type = arg$type,
                                plot_3d = plot_3d[1],
                                plot_sens = TRUE,
                                const = arg$const,
                                compound = arg$compound,### you dont need compund here
                                varlist = sens_varlist,
                                calledfrom =  "plot",
                                npar = arg$npar,
                                calculate_criterion = calculate_criterion,
                                silent = silent,
                                calculate_sens = sensitivity)
  }
  if (evolution){
    ## extract evolution data from the object
    mean_cost <- sapply(1:totaliter, FUN = function(j)obj$evol[[j]]$mean_cost)
    min_cost <- sapply(1:totaliter, FUN = function(j)obj$evol[[j]]$min_cost)
    if (calculate_criterion)
      min_cost[totaliter] <- sens_res$crtval


    type <- obj$arg$type

    # plot setting
    legend_place <- "topright"
    legend_text <- c( "Best Imperialist", "Mean of Imperialists")
    line_col <- c("firebrick3", "blue4")
    title1 <- "Bayesian criterion"
    ylim = switch(type, "bayesian_D" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2))


    PlotEffCost(from = 1,
                to = (totaliter),
                AllCost = min_cost, ##all criterion up to now (all cost function)
                UserCost = NULL,
                DesignType = type,
                col = line_col[1],
                xlab = "Iteration",
                ylim = ylim,
                lty = 1,
                title1 = title1,
                plot_main = TRUE)
    ICAMean_line <-  mean_cost[1:(totaliter)]
    lines(x = 1:(totaliter),
          y = ICAMean_line,
          col = line_col[2], type = "s", lty = 5)
    legend(x = legend_place,  legend = legend_text,lty = c(1,5, 3), col = line_col, xpd = TRUE, bty = "n")
  }
  if(sensitivity || calculate_criterion)
    return(sens_res) else
      return(invisible(NULL))
}
######################################################################################################*
######################################################################################################*
#  roxygen
#' Printing \code{bayes} Objects
#'
#' Print method for an object of class \code{bayes}.
#' @param x an object of class \code{bayes}.
#' @param iter Iteration number. if \code{NULL}, will be set to last iteration.
#' @param ... Argument with no further use.
#' @seealso \code{\link{bayes}}
#' @export

print.bayes <- function(x, iter = NULL, ...){


  if (any(class(x) != c("list", "bayes")))
    stop("'x' must be of class 'bayes'")
  ## to not get confused with design points
  object <- x
  if (is.null(iter))
    totaliter <- length(object$evol) else
      totaliter <- iter
  if (totaliter > length(x$evol))
    stop("'iter' is larger than the maximum number of iterations")
  # if( grepl("on_average", x$arg$type))
  #   type <- "robust" else
  type <- x$arg$type
  ### printing, match with cat in iterate functions
  cat("\n***********************************************************************",
      "\nICA iter:", totaliter, "\n",
      print_xw_char(x = object$evol[[totaliter]]$x, w =  object$evol[[totaliter]]$w, npred = length(object$arg$lx)),
      "\nCriterion value: ", object$evol[[totaliter]]$min_cost,
      "\nTotal number of function evaluations:", object$alg$nfeval,
      "\nTotal number of successful local search moves:", object$alg$nlocal,
      "\nTotal number of successful revolution moves:", object$alg$nrevol,
      "\nConvergence:", object$alg$convergence)
  if (object$arg$ICA.control$only_improve)
    cat("\nTotal number of successful assimilation moves:", object$alg$nimprove, "\n")
  cat("***********************************************************************")
  if (!is.null(object$evol[[totaliter]]$sens))
    print(object$evol[[totaliter]]$sens)
  return(invisible(NULL))
}
######################################################################################################*
######################################################################################################*
#' Printing \code{sensbayes} Objects
#'
#' Print method for an object of class \code{sensbayes}.
#' @param x An object of class \code{sensbayes}.
#' @param ... Argument with no further use.
#' @export
#' @seealso \code{\link{sensbayes}}, \code{\link{sensbayescomp}}
print.sensbayes <- function(x,  ...){
  if (any(class(x) != c("list", "sensbayes")))
    stop("'x' must be of class 'sensbayes'")

  cat("\n***********************************************************************",
      "\nMaximum of the sensitivity function is ", x$max_deriv, "\nEfficiency lower bound (ELB) is ", x$ELB,
      "\nVerification required",x$time, "seconds!",
      "\nAdjust the control parameters in 'sens.bayes.control' for higher speed\n")
  #   if (x$type == "minimax")
  #     optimchar <- "\nSet of all maxima over parameter space\n"
  #   if (x$type == "standardized")
  #     optimchar <- "\nSet of all minima over parameter space\n"
  #   cat("\nAnswering set:\n", answer, optimchar, optima, "\n")
  # }
  #
  # if (!is.null(x$crtval))
  #   cat("Criterion value found by 'nloptr':", x$crtval)
  cat("***********************************************************************")

  return(invisible(NULL))
}
######################################################################################################*
######################################################################################################*
#' @title Control Parameters for Verifying General Equivalence Theorem for Bayesian Designs
#'
#' @description The function \code{sens.bayes.control} returns a list of \code{\link[cubature]{hcubature}}  control parameters for  approximating the integrals in the sensitivity (derivative) function of Bayesian criteria
#' and also \code{\link[nloptr]{nloptr}} control parameters to find maximum of the sensitivity (derivative) function over the design space
#'  and calculate the efficiency lower bound (ELB).
#' @param cubature A list that will be passed to the arguments of the \code{\link[cubature]{hcubature}} function. See 'Details'.
#' @param x0 Vector of starting values for maximizing the sensitivity (derivative) function over the design space \eqn{x}.
#' It will be passed to the optimization function \code{\link[nloptr]{nloptr}}.
#' @param optslist A list will be passed to \code{opts} argument of the function \code{\link[nloptr]{nloptr}} to find the maximum of the sensitivity function over the design space. See 'Details'.
#' @param ... Further arguments will be passed to \code{\link{nl.opts}} from package \code{\link[nloptr]{nloptr}}.
#' @return A list of control parameters for verifying the general equivalence theorem with respect to the Bayesian optimality criteria.
#' @details
#' \code{cubature} is a list that its components will be passed to the function \code{\link[cubature]{hcubature}}.
#' Its components are:
#'  \describe{
#'   \item{\code{tol}}{The maximum tolerance. Defaults to \code{1e-6}.}
#'   \item{\code{maxEval}}{The maximum number of function evaluations needed. Note that the actual number of function evaluations performed is only approximately guaranteed not to exceed this number. Defaults to \code{100000}.}
#'   \item{\code{absError}}{The maximum absolute error tolerated. Defaults to \code{0}.}
#' }
#'
#' ELB is a measure of  proximity of a design to the optimal design without knowing the latter.
#' Given a design, let \eqn{\epsilon} be the global maximum
#'  of the sensitivity (derivative) function with respect the vector of the model predictors \eqn{x} over the design space.
#' ELB is given by \deqn{ELB = p/(p + \epsilon),}
#' where \eqn{p} is the number of model parameters. Obviously,
#' calculating ELB requires finding \eqn{\epsilon} and therefore,
#' a maximization problem to be solved. The function \code{\link[nloptr]{nloptr}}
#' is used here to solve this maximization problem. The arguments \code{x0} and \code{optslist}
#' will be passed to this function as follows:
#'
#' Argument \code{x0} provides the user initial values for this maximization problem
#'  and will be passed to the argument with the same name
#' of the function  \code{\link[nloptr]{nloptr}}.
#'
#'
#' Argument \code{optslist} will be passed to the argument \code{opts} of the function \code{\link[nloptr]{nloptr}}.
#' \code{optslist} is a \code{list} and the most important components are listed as follows:
#'  \describe{
#'   \item{\code{stopval}}{Stop minimization when an objective value <= \code{stopval} is found. Setting stopval to \code{-Inf} disables this stopping criterion (default).}
#'   \item{\code{algorithm}}{Defaults to \code{NLOPT_GN_DIRECT_L}. DIRECT-L is a deterministic-search algorithm based on systematic division of the search domain into smaller and smaller hyperrectangles.}
#'   \item{\code{xtol_rel}}{Stop when an optimization step (or an estimate of the optimum) changes every parameter by less than \code{xtol_rel} multiplied by the absolute value of the parameter. Criterion is disabled if \code{xtol_rel} is non-positive. Defaults to \code{1e-8}.}
#'   \item{\code{ftol_rel}}{Stop when an optimization step (or an estimate of the optimum) changes the objective function value by less than \code{ftol_rel} multiplied by the absolute value of the function value. Criterion is disabled if \code{ftol_rel} is non-positive. Defaults to \code{1e-10}.}
#'   \item{\code{maxeval}}{Stop when the number of function evaluations exceeds maxeval. Criterion is disabled if maxeval is non-positive. Defaults to \code{6000}. See below.}
#' }
#'  A full description of all options is shown by the function \code{nloptr.print.options()} in package \code{\link[nloptr]{nloptr}}.

#' @note  When the value of ELB is larger than 1, it means the maximum found by the optimization function set by \code{algorithm} is not global.
#'  In this case, please increase  the value of the parameter \code{maxeval} to find the global maximum of the sensitivity (derivative) function and avoid false ELB.
#'
#' @export
#' @examples
#' sens.bayes.control()
#' sens.bayes.control(cubature = list(maxEval = 50000))
#' sens.bayes.control(optslist = list(maxeval = 3000))
sens.bayes.control <- function(cubature = list(tol = 1e-6,
                                               maxEval = 100000,
                                               absError = 0),
                               x0 = NULL,
                               optslist = list(stopval = -Inf,
                                               algorithm = "NLOPT_GN_DIRECT_L",
                                               xtol_rel = 1e-8,
                                               ftol_rel = 1e-10,
                                               maxeval = 2000),
                               ...){

  optstlist2 <- do.call(c, list(optslist, list(...)))
  outlist <- suppressWarnings(nloptr::nl.opts(optstlist2))
  outlist["algorithm"] <- optslist$algorithm

  ## outlist has the defaut values of the nl.opts, when any component is null.
  # we play with that here
  if (is.null(optslist$algorithm))
    outlist$algorithm <- "NLOPT_GN_DIRECT_L"
  if (is.null(optslist$stopval))
    outlist$stopval<- -Inf
  if (is.null(optslist$xtol_rel))
    outlist$xtol_rel <- 1e-8
  if (is.null(optslist$ftol_rel))
    outlist$ftol_rel <- 1e-10
  if (is.null(optslist$maxeval))
    outlist$maxeval <- 2000


  ### cubature part
  cubature_out <- do.call(control.cubature, cubature)
  if (is.null(cubature$tol))
    cubature_out$tol <- 1e-6
  if (is.null(cubature$maxEval))
    cubature_out$maxEval <- 100000
  # if (is.null(cubature$doChecking))
  #   cubature_out$doChecking <- FALSE
  # if (is.null(cubature$norm))
  #   cubature_out$norm <- c("INDIVIDUAL", "PAIRED", "L2", "L1", "LINF")
  if (is.null(cubature$absError))
    cubature_out$absError <- 0


  return(list(x0 = x0, optslist = outlist, cubature = cubature_out))
}
######################################################################################################*
######################################################################################################*
#' @title Control Parameters for Evaluating Bayesian Criteria
#'
#' @description The function \code{crt.bayes.control} returns a list of \code{\link[cubature]{hcubature}}  control parameters
#'  for  approximating the integrals in  Bayesian criteria. The key tuning parameters here
#'  are \strong{\code{tol}} and \strong{\code{maxEval}}. Their value affect the algorithm speed and
#'  the accuracy of the results.
#'  The user should find a trade-off between accuracy and speed for his/her example.
#' @param cubature A list that will be passed to the arguments of the function \code{\link[cubature]{hcubature}}. See 'Details'.
#'
#' @details
#' \code{cubature} is a list that its components will be passed to the function \code{\link[cubature]{hcubature}}.
#' Its components are:
#'  \describe{
#'   \item{\code{tol}}{The maximum tolerance. Defaults to \code{1e-5}.}
#'   \item{\code{maxEval}}{The maximum number of function evaluations needed. Note that the actual number of function evaluations performed is only approximately guaranteed not to exceed this number. Defaults to \code{5000}.}
#'   \item{\code{absError}}{The maximum absolute error tolerated. Defaults to \code{0}.}
#' }
#'
#' One can specify a maximum number of function evaluations.
#'  Otherwise, the integration stops when the estimated error is less than
#'   the absolute error requested, or when the estimated error is less than
#'    tol times the integral, in absolute value, or the maximum number of iterations
#'     is reached, whichever is earlier.
#' @examples
#' crt.bayes.control()
#' crt.bayes.control(cubature = list(tol = 1e-4))
#' @return A list of control parameters for \code{\link[cubature]{hcubature}}.
#' @export
crt.bayes.control <- function(cubature = list(tol = 1e-5, maxEval = 50000, absError = 0)){
  cubature_out <- do.call(control.cubature, cubature)
  if (is.null(cubature$tol))
    cubature_out$tol <- 1e-5
  if (is.null(cubature$maxEval))
    cubature_out$maxEval <- 50000
  # if (is.null(cubature$doChecking))
  #   cubature_out$doChecking <- FALSE
  # if (is.null(cubature$norm))
  #   cubature_out$norm <- c("INDIVIDUAL", "PAIRED", "L2", "L1", "LINF")
  if (is.null(cubature$absError))
    cubature_out$absError <- 0
  return(list(cubature = cubature_out))
}
######################################################################################################*
######################################################################################################*

######################################################################################################*
######################################################################################################*
# roxygen
#' @title Updating an Object of Class \code{bayes}
#'
#' @description  Runs the ICA optimization algorithm on an object of class \code{bayes} for more number of iterations  and updates the results.
#'
#' @param object An object of class \code{bayes}.
#' @param iter Number of iterations.
#' @seealso \code{\link{bayes}}
#' @export


# @importFrom nloptr directL you have it in minimax
## @importFrom sn dmsn dmst dmsc
# @importFrom LaplacesDemon dmvl dmvt dmvc dmvpe
iterate.bayes <- function(object, iter){
  if (all(class(object) != c("list", "bayes")))
    stop("''object' must be of class 'bayes'")
  if (missing(iter))
    stop("'iter' is missing")

  arg <- object$arg
  ICA.control <- object$arg$ICA.control
  crt.bayes.control <- object$arg$crt.bayes.control
  sens.bayes.control <-  object$arg$sens.bayes.control
  evol <- object$evol
  npred <- length(arg$lx)
  type <- arg$type
  ## number of parameters
  #npar <- arg$npar
  if (!(type %in% c("D", "DPA", "DPM", "multiple")))
    stop("bug: 'type' must be  'D' or 'DPM' or 'DPM' or 'multiple' in 'iterate.ICAB")
  if (ICA.control$equal_weight)
    w_equal <- rep(1/arg$k, arg$k)

  #############################################################################*
  # plot setting
  #plot_cost <- control$plot_cost
  #plot_sens <- control$plot_sens
  legend_place <- "topright"
  legend_text <- c( "Best Imperialist", "Mean of Imperialists")
  line_col <- c("firebrick3", "blue4")
  title1 <- "Bayesian criterion"

  ################################################################################*

  ## In last iteration the check functions should be applied??
  check_last <- ifelse(ICA.control$checkfreq != FALSE, TRUE, FALSE)

  ################################################################################*
  ### Psi as a function of x and x, y for plotting. Psi_x defined as minus psi to find the minimum
  ## Psi_x is mult-dimensional, x can be of two dimesnion.
  if(length(arg$lx) == 1)
    Psi_x_plot <-  arg$Psi_x ## for PlotPsi_x

  # it is necessary to distniguish between Psi_x for plotiing and finding ELB becasue in plotting for models with two
  # explanatory variables the function should be defined as a function of x, y (x, y here are the ploints to be plotted)

  if(length(arg$lx) == 2)
    Psi_x_plot <- arg$Psi_xy
  #when length(lx) == 1, then Psi_x_plot = Psi_x
  ################################################################################*

  ############################################################################################################*
  ## x_id, w_id are the index of x and w in positions
  #cost_id is the index of
  ## in symmetric case the length of x_id can be one less than the w_id if the number of design points be odd!
  if (ICA.control$sym)
    x_id <- 1:floor(arg$k/2) else
      x_id <- 1:(arg$k * npred)
  if (!ICA.control$equal_weight)
    w_id <- (x_id[length(x_id)] + 1):length(arg$ld) else
      w_id <- NA


  ######################################################################################################*
  ## whenever Calculate_Cost is used, the fixed_arg list should be passed to
  ## fixed argumnet for function Calculate_Cost
  fixed_arg = list(x_id = x_id,
                   w_id = w_id,
                   sym = ICA.control$sym ,
                   sym_point = ICA.control$sym_point,
                   npred = npred,
                   equal_weight = ICA.control$equal_weight,
                   k = arg$k,
                   crfunc = arg$crfunc,
                   Calculate_Cost = Calculate_Cost_bayes)

  vertices_outer <- make_vertices(lower = arg$lx, upper = arg$ux)
  sens_varlist <-list(npred = npred,
                      # plot_3d = "lattice",
                      npar = arg$npar,
                      fimfunc_sens = arg$FIM_sens,
                      Psi_x_bayes  = arg$Psi_funcs$Psi_x_bayes,
                      Psi_xy_bayes  = arg$Psi_funcs$Psi_xy_bayes,
                      crfunc = arg$crfunc,
                      vertices_outer = vertices_outer)
  ########################################################################################*

  #################################################################################################*
  # Initialization when evol is NULL
  #################################################################################################*
  if (is.null(evol)){
    ## set the old seed if call is from minimax
    if (!is.null(ICA.control$rseed))
      set.seed(ICA.control$rseed)
    msg <- NULL
    revol_rate <- ICA.control$revol_rate
    maxiter <- iter
    totaliter <- 0
    #evol <- list()
    min_cost <- c() ## cost of the best imperialists
    mean_cost <- c() ## mean cost of all imperialists
    check_counter <- 0 ## counter to count the check
    total_nlocal  <- 0 ## total number of successful local search
    if (!ICA.control$lsearch)
      total_nlocal <- NA
    total_nrevol <- 0 ## total number of successful revolution
    total_nimprove <- 0 ##total number of improvements due to assimilation

    ############################################## Initialization for ICA
    InitialCountries <- GenerateNewCountry(NumOfCountries = ICA.control$ncount,
                                           lower = arg$ld,
                                           upper = arg$ud,
                                           sym = ICA.control$sym,
                                           w_id = w_id,
                                           x_id = x_id,
                                           npred= npred,
                                           equal_weight = ICA.control$equal_weight)
    if (!is.null(arg$initial))
      InitialCountries[1:dim(arg$initial)[1], ] <- arg$initial
    InitialCost <- vector("double", ICA.control$ncount)

    temp <- fixed_arg$Calculate_Cost(mat = InitialCountries, fixed_arg = fixed_arg)

    total_nfeval <-  temp$nfeval
    InitialCost <-  temp$cost
    inparam <- temp$inner_optima ## we require that to avoid errors!!
    ## waring inparam for optim_on_average does not have any meaning!
    temp <- NA # safety
    ##Now we should sort the initial countries with respect to their initial cost
    SortInd <- order(InitialCost)
    InitialCost <- InitialCost[SortInd] # Sort the cost in assending order. The best countries will be in higher places
    InitialCountries <- InitialCountries[SortInd,, drop = FALSE] #  Sort the population with respect to their cost. The best country is in the first column


    # creating empires
    Empires <- CreateInitialEmpires(sorted_Countries = InitialCountries,
                                    sorted_Cost = InitialCost,
                                    Zeta = ICA.control$zeta,
                                    sorted_InnerParam = inparam,
                                    NumOfInitialImperialists = ICA.control$nimp,
                                    NumOfAllColonies = (ICA.control$ncount -ICA.control$nimp))

    best_imp_id<- 1 ## the index of list in which best imperialists is in.

    ########################################################################*
  }
  #################################################################################################*

  #################################################################################################*
  # when we are updating the object for more number of iterations
  #################################################################################################*
  if (!is.null(evol)){
    ## reset the seed!
    # if (exists(".Random.seed")){
    #   GlobalSeed <- get(".Random.seed", envir = .GlobalEnv)
    #   #if you call directly from iterate and not minimax!
    #   on.exit(assign(".Random.seed", GlobalSeed, envir = .GlobalEnv))
    # }
    msg <- object$best$msg
    prev_iter <- length(evol) ##previous number of iterationst
    maxiter <- iter + prev_iter
    totaliter <- prev_iter
    mean_cost <- sapply(1:(totaliter), FUN = function(j) evol[[j]]$mean_cost)
    min_cost <- sapply(1:(totaliter), FUN = function(j) evol[[j]]$min_cost)
    Empires <- object$empires


    check_counter <- arg$updating$check_counter
    total_nfeval <- object$alg$nfeval
    total_nlocal <-  object$alg$nlocal
    total_nrevol <-object$alg$nrevol
    total_nimprove <-object$alg$nimprove

    imp_cost <- round(sapply(object$empires, "[[", "ImperialistCost"), 12)
    best_imp_id<- which.min(imp_cost)
    revol_rate <-  arg$updating$revol_rate
    ##updating the random seed
    # if (!is.null(ICA.control$rseed)){
    #   do.call("RNGkind",as.list(arg$updating$oldRNGkind))  ## must be first!
    #   assign(".Random.seed", arg$updating$oldseed , .GlobalEnv)
    # }
  }
  ##########################################################################*
  space_size <- arg$ud - arg$ld
  continue = TRUE
  vertices_outer <- make_vertices(lower = arg$lx, upper = arg$ux) ## we need it for checking the equivalence theorem
  #################################################################################################################*
  ### start of the while loop until continue == TRUE
  #################################################################################################################*
  while (continue == TRUE){
    totaliter <- totaliter + 1


    check_counter <- check_counter + 1
    revol_rate <- ICA.control$damp * revol_rate
    ## revolution rate is increased by damp ration in every iter


    ###############################################################################################*
    ################################################################# for loop over all empires[ii]
    for(ii in 1:length(Empires)){
      #cat(totaliter, " while loop: ", ii, "\n")
      ########################################## local search is only for point!
      if (ICA.control$lsearch){
        LocalSearch_res <- LocalSearch (TheEmpire =  Empires[[ii]],
                                        lower = arg$ld,
                                        upper = arg$ud,
                                        l = ICA.control$l,
                                        fixed_arg = fixed_arg)

        Empires[[ii]] <- LocalSearch_res$TheEmpire
        total_nfeval <- total_nfeval + LocalSearch_res$nfeval
        total_nlocal <- total_nlocal + LocalSearch_res$n_success
      }
      ##########################################################################*
      # if (totaliter == 2 & ii == 4)
      #   debug(Calculate_Cost_bayes)
      ############################################################## Assimilation
      temp5 <- AssimilateColonies2(TheEmpire = Empires[[ii]],
                                   AssimilationCoefficient = ICA.control$assim_coeff,
                                   VarMin = arg$ld,
                                   VarMax = arg$ud,
                                   ExceedStrategy = "perturbed",
                                   sym = ICA.control$sym,
                                   AsssimilationStrategy = ICA.control$assim_strategy,
                                   MoveOnlyWhenImprove = ICA.control$only_improve,
                                   fixed_arg = fixed_arg,
                                   w_id = w_id,
                                   equal_weight = ICA.control$equal_weight)
      ##Warning: in this function the colonies position are changed but the imperialist and the
      ##cost functions of colonies are not updated yet!
      ##they will be updated after revolution
      Empires[[ii]] <- temp5$TheEmpire
      total_nfeval <- total_nfeval + temp5$nfeval
      total_nimprove <-  total_nimprove + temp5$nimprove
      ##########################################################################*

      ############################################################### Revolution
      temp4 <- RevolveColonies(TheEmpire = Empires[[ii]],
                               RevolutionRate = revol_rate,
                               NumOfCountries = ICA.control$ncount,
                               lower = arg$ld,
                               upper = arg$ud,
                               sym = ICA.control$sym,
                               sym_point = ICA.control$sym_point,
                               fixed_arg = fixed_arg,
                               w_id = w_id,
                               equal_weight = ICA.control$equal_weight)
      Empires[[ii]] <- temp4$TheEmpire
      total_nrevol <- total_nrevol + temp4$nrevol
      total_nfeval <- total_nfeval + temp4$nfeval
      ############################################################*
      Empires[[ii]] <- PossesEmpire(TheEmpire = Empires[[ii]])

      ##after updating the empire the total cost should be updated
      ## Computation of Total Cost for Empires
      Empires[[ii]]$TotalCost <- Empires[[ii]]$ImperialistCost + ICA.control$zeta * mean(Empires[[ii]]$ColoniesCost)

    }
    ############################################################ end of the loop for empires [[ii]]
    ###############################################################################################*

    #################################################### Uniting Similiar Empires
    if (length(Empires)>1){

      Empires <- UniteSimilarEmpires(Empires = Empires,
                                     Zeta = ICA.control$zeta,
                                     UnitingThreshold = ICA.control$uniting_threshold,
                                     SearchSpaceSize = space_size)
    }
    ############################################################################*
    # zeta is necessary to update the total cost!
    Empires <- ImperialisticCompetition(Empires = Empires, Zeta = ICA.control$zeta)

    ############################################################## save the seed
    # we get the seed here because we dont know if in cheking it wil be chaned
    #we save the seed when we exit the algorithm
    oldseed <- get(".Random.seed", envir = .GlobalEnv)
    oldRNGkind <- RNGkind()
    ############################################################################*

    ############################################################################*
    # extracing the best emperor and its position
    imp_cost <- round(sapply(Empires, "[[", "ImperialistCost"), 12)
    if (type %in% c("D", "DPA", "DPM", "multiple")){
      min_cost[totaliter] <-   min(imp_cost)
      mean_cost[totaliter] <-  mean(imp_cost)
    }else
      stop("Bug: check the type")
    best_imp_id <- which.min(imp_cost) ## which list contain the best imp
    if (!ICA.control$equal_weight)
      w <- Empires[[best_imp_id]]$ImperialistPosition[, w_id] else
        w <- w_equal
    x <- Empires[[best_imp_id]]$ImperialistPosition[, x_id]



    if (ICA.control$sym){
      x_w <- ICA_extract_x_w(x = x, w = w, sym_point = ICA.control$sym_point)
      x <- x_w$x
      w <- x_w$w
    }
    ##sort Point
    if (npred == 1){
      w <- w[order(x)]
      x <- sort(x)
    }
    ############################################################################*

    ################################################################ print trace

    if (ICA.control$trace){
      cat("\nIteration:", totaliter, "\nDesign Points:\n", x, "\nWeights: \n", w,
          "\nCriterion value: ", min_cost[totaliter],
          "\nTotal number of function evaluations:", total_nfeval, "\nTotal number of successful local search moves:", total_nlocal,
          "\nTotal number of successful revolution moves:", total_nrevol, "\n")
      if (ICA.control$only_improve)
        cat("Total number of successful assimilation moves:", total_nimprove, "\n")

    }
    ############################################################################*
    if ( min_cost[totaliter] == 1e-24)
      warning("Computational issue! maybe the design is singular!\n")

    ################################################################### continue
    if (totaliter ==  maxiter){
      continue <- FALSE
      convergence = "maxiter"
    }
    if(length(Empires) == 1 && ICA.control$stop_rule == "one_empire"){
      continue <- FALSE
      convergence = "one_empire"
    }
    ## the continue also can be changed in check
    ############################################################################*

    ################################################################################# plot_cost
    if (ICA.control$plot_cost) {
      #ylim for efficiency depends on the criterion type because
      ylim = switch(type, "D" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2))
      PlotEffCost(from = 1,
                  to = (totaliter),
                  AllCost = min_cost, ##all criterion up to now (all cost function)
                  UserCost = NULL,
                  DesignType = type,
                  col = line_col[1],
                  xlab = "Iteration",
                  ylim = ylim,
                  lty = 1,
                  title1 = title1,
                  plot_main = TRUE)
      ICAMean_line <-  mean_cost[1:(totaliter)]
      lines(x = 1:(totaliter),
            y = ICAMean_line,
            col = line_col[2], type = "s", lty = 5)
      legend(x = legend_place,  legend = legend_text,lty = c(1,5, 3), col = line_col, xpd = TRUE, bty = "n")
    }
    ############################################################################*

    ####################################################################################*
    #check the equivalence theorem and find ELB
    ####################################################################################*
    ## we check the quvalence theorem in the last iteration anyway. but we may not plot it.
    if (check_counter == ICA.control$checkfreq || (check_last && !continue)){
      if (arg$ICA.control$trace){
        #cat("\n*********************************************************************")
        if (!continue)
          cat("\nOptimization is done!\n")
        cat("Requesting design verification by the general equivalence theorem\n")
      }
      sens_res <- sensbayes_inner(x = x, w = w, lx = arg$lx, ux = arg$ux,
                                  fimfunc = arg$FIM, prior = arg$prior,
                                  sens.bayes.control = sens.bayes.control,
                                  crt.bayes.control = crt.bayes.control,
                                  type = arg$type,
                                  plot_3d = arg$plot_3d,
                                  plot_sens = ICA.control$plot_sens,
                                  const = arg$const,
                                  compound = arg$compound,
                                  varlist = sens_varlist,
                                  calledfrom = "iter",
                                  npar = arg$npar,
                                  calculate_criterion = FALSE,
                                  silent = !arg$ICA.control$trace)

      GE_confirmation <- (sens_res$ELB >= ICA.control$stoptol)
      ##########################################################################*
      # print trace that is related to checking
      # if (ICA.control$trace)
      #   cat("maximum of sensitivity:", sens_res$max_deriv, "\nefficiency lower bound (ELB):", sens_res$ELB, "\n")
      ##########################################################################*

      #if (npred == 1){
      if (GE_confirmation && ICA.control$stop_rule == "equivalence"){
        continue <- FALSE
        convergence <- "equivalence"
      }
      ##########################################################################*
    }else
      sens_res <- NULL
    ####################################################################### end of check
    #  if (check_counter == control$checkfreq || (check_last && !continue)) ##########
    ####################################################################################*

    ####################################################################### save
    evol[[totaliter]] <- list(iter = totaliter, x = x, w = w, min_cost = min_cost[totaliter], mean_cost = mean_cost[totaliter], sens = sens_res)
    ############################################################################*

    ################################################################ print trace
    # if (ICA.control$trace){
    #   cat("total local search:", total_nlocal, "\n")
    #   cat("total revolution:", total_nrevol, "\n")
    #   if (ICA.control$only_improve)
    #     cat("total improve:", total_nimprove, "\n")
    # }
    ############################################################################*
  }
  #################################################################################################################*
  ### end of the while loop over continue == TRUE
  #################################################################################################################*

  if (!ICA.control$only_improve)
    total_nimprove <- NA
  msg <- NULL
  ##############################################################################*

  ######################################################################## saving
  ## we add the following to arg becasue dont want to document it in Rd files
  # updating parameters
  object$arg$updating$check_counter <- check_counter
  object$arg$updating$oldseed <- oldseed
  object$arg$updating$oldRNGkind <-  oldRNGkind
  object$arg$updating$revol_rate = revol_rate ## different from revolrate

  object$evol <- evol
  object$empires <- Empires
  object$alg <- list(
    nfeval = total_nfeval,
    nlocal = total_nlocal,
    nrevol = total_nrevol,
    nimprove = total_nimprove,
    convergence = convergence)

  ## arg has a list named update as well
  ###### end of saving
  ##############################################################################*
  return(object)

}


######################################################################################################*
######################################################################################################*



######################################################################################################*
######################################################################################################*
