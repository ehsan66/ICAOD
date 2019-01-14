######################################################################################################*
######################################################################################################*
#' @title  Bayesian D-Optimal Designs
#'
#' @description
#'  Finds (pseudo) Bayesian D-optimal designs for nonlinear models.
#'  It should be used when the user assumes a (truncated) prior distribution for the unknown model parameters.
#'  If you have a discrete prior, please use the function \code{\link{robust}}.
#'
#' @inheritParams minimax
#' @param prior An object of class \code{cprior}. User can also use one of the functions
#'  \code{\link{uniform}}, \code{\link{normal}},
#' \code{\link{skewnormal}} or \code{\link{student}}  to create the  prior. See 'Details' of \code{\link{bayes}}.
#' @param crt.bayes.control A list. Control parameters to approximate the integral in  the Bayesian criterion at a given design over the parameter space.
#'  For details, see \code{\link{crt.bayes.control}}.
#' @param sens.bayes.control A list. Control parameters to verify the general equivalence theorem. For details, see \code{\link{sens.bayes.control}}.
#' @param crt_method A character denotes which method to be used to approximate the integrals in Bayesian criteria.
#'  \code{"cubature"} corresponds to the adaptive multivariate integration method using the \code{\link[cubature]{hcubature}} algorithm (default).
#'  \code{"quadrature"} corresponds the traditional quadrature formulas and calls the function \code{\link[mvQuad]{createNIGrid}}.
#'  The tuning parameters are adjusted by \code{crt.bayes.control}. Default is set to \code{"cubature"}.
#' @param sens_method Similar to \code{crt_method}, but for approximating the integrals in the sensitivity (derivative) function.
#' The tuning parameters are adjusted by \code{sens.bayes.control}. Default is set to \code{"cubature"}.
#' @param npar Number of model parameters.  Used when \code{fimfunc} is given instead of \code{formula} to specify the number of model parameters.
#'   If not specified correctly, the sensitivity (derivative) plot may be shifted below the y-axis.
#'   When \code{NULL} (default), it will be set to \code{length(parvars)} or
#'   \code{prior$npar} when \code{missing(formula)}.
#' @export
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
#'  \item{\code{fn}: }{Prior distribution as an R \code{function} with argument \code{param} that is the vector of the unknown parameters. See below.}
#'  \item{\code{npar}: }{Number of unknown parameters and is equal to the length of \code{param}}.
#'  \item{\code{lower}: }{Argument \code{lower}. It has the same length as \code{param}}.
#'  \item{\code{upper}: }{Argument \code{upper}. It has the same length as \code{param}}.
#' }
#' A \code{cprior} object  will be passed to the argument \code{prior} of the function \code{\link{bayes}}.
#'  The argument \code{param} in \code{fn} has the same order as the argument \code{parvars} when the model is specified by a \code{formula}.
#' Otherwise, it is the same as the argument \code{param} in the function \code{fimfunc}.\cr
#' The user can use the implemented  priors that are \code{\link{uniform}}, \code{\link{normal}},
#' \code{\link{skewnormal}} and \code{\link{student}} to create a \code{cprior} object.
#'
#' To verify the equivalence theorem of the output design,
#' use \code{\link{plot}} function or change the value of the \code{checkfreq} in the argument \code{\link{ICA.control}}.
#'
#' To increase the speed of the algorithm, change the value of the tuning parameters \code{tol} and \code{maxEval} via
#' the argument  \code{crt.bayes.control} when \code{crt_method = "cubature"}.
#' Similarly, this applies  when \code{crt_method = "quadrature"}.
#' In general, if the CPU time matters, the user should find an appropriate speed-accuracy trade-off  for her/his own problem.
#'  See 'Examples' for more details.
#'
#' If some of the parameters are known and fixed, they should be set
#' to their values via the argument \code{paravars} when the model is given by \code{formula}. In this case,
#' the user must provide the number of parameters via the argument \code{npar} for verifying the general equivalence theorem.
#'  See 'Examples', Section 'Weibull',  'Richards' and 'Exponential' model.
#' @return
#'  an object of class \code{bayes} that is a list including three sub-lists:
#' \describe{
#'   \item{\code{arg}}{A list of design and algorithm parameters.}
#'   \item{\code{evol}}{A list of length equal to the number of iterations that stores the information about the best design (design with the minimum criterion value) of each iteration as follows:
#'    \code{evol[[iter]]} contains:
#'     \tabular{lll}{
#'       \code{iter}                   \tab      \tab Iteration number.\cr
#'       \code{x}                      \tab      \tab Design points. \cr
#'       \code{w}                      \tab      \tab Design weights. \cr
#'       \code{min_cost}               \tab      \tab Value of the criterion for the best imperialist (design).  \cr
#'       \code{mean_cost}              \tab      \tab Mean of the criterion values of all the imperialists. \cr
#'       \code{sens}                   \tab      \tab An object of class \code{'sensbayes'}. See below.\cr
#'     }
#'   }
#'
#'   \item{\code{empires}}{A list of all the empires of the last iteration.}
#'   \item{\code{alg}}{A list with following information:
#'     \tabular{lll}{
#'       \code{nfeval}           \tab      \tab Number of function evaluations. It does not count the function evaluations from checking the general equivalence theorem. \cr
#'       \code{nlocal}           \tab      \tab Number of successful local searches. \cr
#'       \code{nrevol}           \tab      \tab Number of successful revolutions. \cr
#'       \code{nimprove}         \tab      \tab Number of successful movements toward the imperialists in the assimilation step. \cr
#'       \code{convergence}      \tab      \tab Stopped by \code{'maxiter'} or \code{'equivalence'}?\cr
#'     }
#'   }
#' }
#' The list \code{sens}  contains information about the design verification by the general equivalence theorem.
#'  See \code{sensbayes} for more Details.
#'   It is only given every \code{ICA.control$checkfreq} iterations
#'  and also the last iteration if   \code{ICA.control$checkfreq >= 0}. Otherwise, \code{NULL}.
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
                  crt_method = c("cubature", "quadrature"),
                  sens_method = c("cubature", "quadrature"),
                  initial = NULL,
                  npar = NULL,
                  plot_3d = c("lattice", "rgl"),
                  x = NULL) {
  #cat("bayes:", get(".Random.seed")[2], "\n")
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
                         crt_method = crt_method[1],
                         sens_method = sens_method[1],
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
                         const = list(ui = NULL, ci = NULL, coef = NULL),
                         only_w_varlist = list(x = x))
  return(output)
}
######################################################################################################*
######################################################################################################*
#' @title Verifying Optimality of Bayesian D-optimal Designs
#' @inheritParams bayes
#' @inheritParams sensminimax
#' @description
#'  Plots the sensitivity (derivative) function  and calculates the efficiency lower bound (ELB) for a  given  Bayesian design.
#' Let \eqn{\boldsymbol{x}}{x} belongs to \eqn{\chi} that denotes the design space.
#' Based on the general equivalence theorem,  a design \eqn{\xi^*}{\xi*} is optimal if and only if the value of the sensitivity (derivative) function
#' is non-positive for all \eqn{\boldsymbol{x}}{x} in \eqn{\chi} and zero when
#'  \eqn{\boldsymbol{x}}{x} belongs to the support of \eqn{\xi^*}{\xi*} (be equal to the one of the design points).
#'
#' For an approximate (continuous) design, when the design space is one or two-dimensional, the user can visually verify the optimality of the design by observing the
#' sensitivity plot. Furthermore, the approximity of the design to the optimal design
#'  can be measured by the  ELB without knowing the latter.
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
#' Depending on the complexity of the problem at hand, sometimes, the CPU time can be considerably reduced
#' by choosing a set of  less conservative values for the tuning parameters \code{tol} and \code{maxEval} in
#' the function \code{\link{sens.bayes.control}} when \code{sens_method = "cubature"}.
#' Similarly, this applies  when \code{sens_method = "quadrature"}.
#' In general, if the CPU time matters, the user should find an appropriate speed-accuracy trade-off  for her/his own problem.
#'  See 'Examples' for more details.
#'
#' @note
#' The default values of the tuning parameters in \code{sens.bayes.control} are set in a way that
#' having accurate plots for the sensitivity (derivative) function
#'  and calculating the ELB to a high precision  to be the primary goals,
#'  although the process may take too long. The user should choose a set of less conservative values
#'  via the argument \code{sens.bayes.control} when the CPU-time is too long or matters.
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
                      crt_method = c("cubature", "quadrature"),
                      sens_method = c("cubature", "quadrature"),
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
                             crt_method = crt_method[1],
                             sens_method = sens_method[1],
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
#'  Finds compound Bayesian DP-optimal designs that meet the dual goal of parameter estimation and
#'   increasing the probability of a particular outcome in a binary response  model.
#'A compound Bayesian DP-optimal design maximizes  the product of the Bayesian efficiencies of a design \eqn{\xi} with respect to D- and average P-optimality, weighted by a pre-defined mixing constant
#' \eqn{0 \leq \alpha \leq 1}{0 <= \alpha <= 1}.
#'
#' @inheritParams bayes
#' @param alpha A value between 0 and 1.
#' Compound or combined DP-criterion  is the product of the efficiencies of a design  with respect to D- and average P- optimality, weighted by \code{alpha}.
#' @param prob Either \code{formula} or a \code{function}. When function, its argument are \code{x} and \code{param}, and they are the same as the arguments in \code{fimfunc}.
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
#'   A  compound Bayesian DP-optimal design maximizes over \eqn{\Xi}
#' \deqn{\int_{\theta \in \Theta} \frac{\alpha}{q}\log|M(\xi, \theta)| + (1- \alpha)
#'\log \left( \sum_{i=1}^k w_ip(x_i, \theta) \right) \pi(\theta) d\theta.}{
#' integration over \Theta \alpha/q log|M(\xi, \theta)| + (1- \alpha)
#'log ( \sum w_i p(x_i, \theta)) \pi(\theta) d\theta.
#'}
#'
#' To verify the equivalence theorem of the output design,
#' use \code{\link{plot}} function or change the value of the \code{checkfreq} in the argument \code{\link{ICA.control}}.
#'
#' To increase the speed of the algorithm, change the value of the tuning parameters \code{tol} and \code{maxEval} via
#' the argument  \code{crt.bayes.control} when \code{crt_method = "cubature"}.
#' In general, if the CPU time matters, the user should find an appropriate speed-accuracy trade-off  for her/his own problem.
#'  See 'Examples' for more details.
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
                      crt_method = c("cubature", "quadrature"),
                      sens_method = c("cubature", "quadrature"),
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
    if (!all(c("x", "param") %in% formalArgs(prob)))
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
                         crt_method = crt_method[1],
                         sens_method = sens_method[1],
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
#'  This function plot the sensitivity (derivative) function given an approximate (continuous) design and calculate the efficiency lower bound (ELB) for Bayesian DP-optimal designs.
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
#' Depending on the complexity of the problem at hand, sometimes, the CPU time can be considerably reduced
#' by choosing a set of  less conservative values for the tuning parameters \code{tol} and \code{maxEval} in
#' the function \code{\link{sens.bayes.control}} when \code{sens_method = "cubature"}. Similarly, this applies  when \code{sens_method = "quadrature"}.
#' In general, if the CPU time matters, the user should find an appropriate speed-accuracy trade-off  for her/his own problem.
#'  See 'Examples' for more details.
#'
#' @note
#' The default values of the tuning parameters in \code{sens.bayes.control} are set in a way that
#' having accurate plots for the sensitivity (derivative) function
#'  and calculating the ELB to a high precision  to be the primary goals,
#'  although the process may take too long. The user should choose a set of less conservative values
#'  via the argument \code{sens.bayes.control} when the CPU-time is too long or matters.
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
                          crt_method = c("cubature", "quadrature"),
                          sens_method = c("cubature", "quadrature"),
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
    if (!all(c("x", "param") %in% formalArgs(prob)))
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
                             crt_method = crt_method[1],
                             sens_method = sens_method[1],
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

######################################################################################################*
######################################################################################################*

#  roxygen
#' Plotting \code{bayes} Objects
#' @description
#'  This function plots the evolution of the ICA algorithm (iteration vs the best (minimum) criterion value at each iteration) and also verifies the optimality of the last obtained design
#'  using  the general equivalence theorem. It plots the sensitivity function and calculates the ELB for the best design generated  at iteration number \code{iter}.
#'
#' @param x An object of class \code{bayes}.
#' @param iter Iteration number. if \code{NULL} (default), it will be set to the last iteration.
#' @param sensitivity Logical. If \code{TRUE} (default), the general equivalence theorem is used to check the optimality if the best design in iteration number \code{iter} and the sensitivity function will be plotted.
#' @param calculate_criterion Logical. Re-calculate the Bayesian criterion value (maybe with a set of new tuning parameters to be sure of the accuracy of the integral approximation)? Defaults to \code{FALSE}. See 'Details'.
#' @param sens.bayes.control Control parameters to verify general equivalence theorem. For details, see \code{\link{sens.bayes.control}}.
#'  If \code{NULL} (default), it will be set to the  tuning parameters used to create object \code{x}.
#' @param crt.bayes.control Control parameters to approximate the integration in the Bayesian criterion at a given design.
#'  For details, see \code{\link{crt.bayes.control}}.
#' If \code{NULL} (default), it will be set to the  tuning parameters used to create object \code{x}.
#' @param silent Do not print anything? Defaults to \code{FALSE}.
#' @param plot_3d Which package should be used to plot the sensitivity function for two-dimensional design space. Defaults to \code{plot_3d = "lattice"}.
#' Only applicable when \code{sensitivity = TRUE}.
#' @param evolution Plot Evolution? Defaults to \code{FALSE}.
#' @param crt_method A character denotes which method to be used to approximate the integrals in Bayesian criteria.
#'  \code{"cubature"} corresponds to the adaptive multivariate integration method using the \code{\link[cubature]{hcubature}} algorithm (default).
#'  \code{"quadrature"} corresponds the traditional quadrature formulas and calls the function \code{\link[mvQuad]{createNIGrid}}.
#'  The tuning parameters are adjusted by \code{crt.bayes.control}.
#'   If \code{NULL} (default), it will use the method already applied to create the object \code{x}.
#' @param sens_method Similar to \code{crt_method}, but for approximating the integrals in the sensitivity (derivative) function.
#' The tuning parameters are adjusted by \code{sens.bayes.control}.
#'  If \code{NULL} (default), it will use the method already applied to create the object \code{x}.
#' @param ... Argument with no further use.
#' @seealso \code{\link{bayes}}, \code{\link{bayescomp}}
#' @details
#' In addition to verifying the general equivalence theorem,
#'  this function makes it possible to  re-calculate the criterion value for the output designs using a new (say, more conservative) set of tuning parameters.
#'  This is useful for  Bayesian optimal designs to assess the robustness of the
#'  criterion value with respect to different values of the tuning parameters.
#'  To put it simple, for Bayesian generated designs, the user can re-calculate the
#'  criterion value for the output design (best design generated in \code{iter}) with different values for  \code{maxEval} and \code{tol} in \code{\link{crt.bayes.control}}
#'  to verify  that the function \code{hcubature} approximates the integrals to an user-acceptable accuracy. The same also applies for the quadrature methods using different number of nodes.
#' @export
plot.bayes <- function(x, iter = NULL,
                       sensitivity = TRUE,
                       calculate_criterion = FALSE,
                       crt_method = NULL,
                       sens_method = NULL,
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

  #### quadrature and cubature
  if (is.null(sens_method))
    sens_method  <- arg$sens_method
  if (is.null(crt_method))
    crt_method  <- arg$crt_method
  if (!crt_method %in% c("cubature", "quadrature"))
    stop("'crt_method' can only be 'cubature' or 'quadrature'")
  if (!sens_method %in% c("cubature", "quadrature"))
    stop("'sens_method' can only be 'cubature' or 'quadrature'")
  ####


  if (calculate_criterion || sensitivity){
    if (is.null(sens.bayes.control)){
      sens.bayes.control <- arg$sens.bayes.control }else {
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
                                   compound = arg$compound,
                                   method =  sens_method[1])
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
                                sens_method = sens_method,
                                crt_method = crt_method,
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
#' @param iter Iteration number. if \code{NULL}, it will be set to the last iteration.
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
      print_xw_char(x = object$evol[[totaliter]]$x, w =  object$evol[[totaliter]]$w, npred = length(object$arg$lx), is.only.w = object$arg$is.only.w),
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
#' @description
#' This function returns two lists each corresponds
#'  to an implemented integration method for approximating the integrals
#'   in the sensitivity (derivative) functions for Bayesian criteria.
#' Moreover, it contains a list of \code{\link[nloptr]{nloptr}} control parameters to find maximum of the sensitivity (derivative) function over the design space,
#'  used to calculate the efficiency lower bound (ELB).
#' @param cubature A list that will be passed to the arguments of the \code{\link[cubature]{hcubature}} function. See 'Details' of \code{\link{crt.bayes.control}}.
#' @param quadrature A list that will be passed to the arguments of the \code{\link[mvQuad]{createNIGrid}} function. See 'Details' of \code{\link{crt.bayes.control}}.
#' @param x0 Vector of starting values for maximizing the sensitivity (derivative) function over the design space \eqn{x}.
#' It will be passed to the optimization function \code{\link[nloptr]{nloptr}}.
#' @param optslist A list will be passed to \code{opts} argument of the function \code{\link[nloptr]{nloptr}} to find the maximum of the sensitivity function over the design space. See 'Details'.
#' @param ... Further arguments will be passed to \code{\link{nl.opts}} from package \code{\link[nloptr]{nloptr}}.
#' @return A list of three lists each contains the  control parameters for verifying the general equivalence theorem with respect to the Bayesian optimality criteria.
#' @details
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
#' in the function  \code{\link[nloptr]{nloptr}}.
#'
#'
#' Argument \code{optslist} is a list and it will be passed to the argument \code{opts} of the function \code{\link[nloptr]{nloptr}}.
#' The most important components of  \code{optslist} are:
#'  \describe{
#'   \item{\code{stopval}}{Stop minimization when an objective value <= \code{stopval} is found. Setting \code{stopval} to \code{-Inf} disables this stopping criterion (default).}
#'   \item{\code{algorithm}}{Defaults to \code{NLOPT_GN_DIRECT_L}. DIRECT-L is a deterministic-search algorithm based on systematic division of the search domain into smaller and smaller hyperrectangles.}
#'   \item{\code{xtol_rel}}{Stop when an optimization step (or an estimate of the optimum) changes every parameter by less than \code{xtol_rel} multiplied by the absolute value of the parameter. Criterion is disabled if \code{xtol_rel} is non-positive. Defaults to \code{1e-8}.}
#'   \item{\code{ftol_rel}}{Stop when an optimization step (or an estimate of the optimum) changes the objective function value by less than \code{ftol_rel} multiplied by the absolute value of the function value. Criterion is disabled if \code{ftol_rel} is non-positive. Defaults to \code{1e-10}.}
#'   \item{\code{maxeval}}{Stop when the number of function evaluations exceeds \code{maxeval}. Criterion is disabled if \code{maxeval} is non-positive. Defaults to \code{6000}. See 'Note' on when to change its value.}
#' }
#'  More details are available by the function \code{nloptr.print.options()} in package \code{\link[nloptr]{nloptr}}.

#' @note
#'  Whenever the value of ELB is computationally larger than 1, it is a sign that  the obtained \eqn{\epsilon} is not global maximum.
#'  To overcome this issue, please increase  the value of the parameter \code{maxeval} to allow the optimization algorithm to find the global maximum of the sensitivity (derivative) function over the design space.
#'
#'
#' Please note the difference between \code{maxEval} in \code{cubature} and \code{maxeval} in  \code{optslist}. They are quite two types of different tuning parameters.
#' The earlier is the (approximate) maximum number of function evaluations in the hcubature adaptive integration method and the latter is the maximum number of function evaluations in the maximization problem.
#' @export
#' @examples
#' sens.bayes.control()
#' sens.bayes.control(cubature = list(maxEval = 50000))
#' sens.bayes.control(optslist = list(maxeval = 3000))
#' sens.bayes.control(quadrature =  list(level = 4))
sens.bayes.control <- function(cubature = list(tol = 1e-5,
                                               maxEval = 50000,
                                               absError = 0),
                               quadrature = list(type = NULL,
                                                 level = NULL,
                                                 ndConstruction = "product",
                                                 level.trans = FALSE),
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
  if (is.null(cubature$absError))
    cubature_out$absError <- 0

  quadrature_out <- do.call(control.quadrature, quadrature)
  if (is.null(quadrature$type))
    quadrature_out$type <- "GLe"
  if (is.null(quadrature$level))
    quadrature_out$level <- 8
  if (is.null(quadrature$lndConstruction))
    quadrature_out$ndConstruction <- "product"
  if (is.null(quadrature$level.trans))
    quadrature_out$level.trans <- FALSE

  return(list(x0 = x0, optslist = outlist, cubature = cubature_out, quadrature = quadrature_out))
}
######################################################################################################*
######################################################################################################*
#' @title Control Parameters for Approximating Bayesian Criteria
#'
#' @description
#'  This function returns two lists each corresponds
#'  to an implemented integration method for approximating the integrals
#'   in Bayesian criteria.
#'  The first list is named \code{cubature} and contains the \code{\link[cubature]{hcubature}}
#'   control parameters to  approximate the integrals with an adaptive multivariate integration method over hypercubes.
#'  The second list is named \code{quadrature} and contains the \code{\link[mvQuad]{createNIGrid}}
#'  tuning parameters to approximate the integrals with the quadrature methods.
#' @param cubature A list that will be passed to the arguments of the function \code{\link[cubature]{hcubature}} for the adaptive multivariate integration.
#'  It is required and used when \code{crt_method = "cubature"} in the parent function, e.g.  \code{\link{bayes}}. See 'Details'.
#' @param quadrature A list that will be passed to the arguments of the function \code{\link[mvQuad]{createNIGrid}} for the quadrature-based integration.
#' It is required and used when \code{crt_method = "quadrature"} in the parent function, e.g.  \code{\link{bayes}}. See 'Details'.
#' @details
#'
#' \code{cubature} is a list that its components will be passed to the function \code{\link[cubature]{hcubature}} and they are:
#'  \describe{
#'   \item{\code{tol}}{The maximum tolerance. Defaults to \code{1e-5}.}
#'   \item{\code{maxEval}}{The maximum number of function evaluations needed. Note that the actual number of function evaluations performed is only approximately guaranteed not to exceed this number. Defaults to \code{5000}.}
#'   \item{\code{absError}}{The maximum absolute error tolerated. Defaults to \code{0}.}
#' }
#'
#' One can specify a maximum number of function evaluations.
#'  Otherwise, the integration stops when the estimated error is less than
#'   the absolute error requested, or when the estimated error is less than
#'    \code{tol} times the absolute value of the integral,  or when the maximum number of iterations
#'     is reached, whichever is earlier.
#'      \code{cubature} is activated when \code{crt_method = "cubature"} in
#'       any of the parent functions (for example, \code{\link{bayes}}).
#'
#' \code{quadrature} is a list that its components will be passed to
#'  the function \code{\link[mvQuad]{createNIGrid}} and they are:
#'  \describe{
#'   \item{\code{type}}{Quadrature rule (see Details of \code{\link[mvQuad]{createNIGrid}}) Defaults to \code{"GLe"}.}
#'   \item{\code{level}}{Accuracy level (typically number of grid points for the underlying 1D quadrature rule). Defaults to \code{6}.}
#'   \item{\code{ndConstruction}}{Character vector which denotes the construction rule
#'    for multidimensional grids. \code{"product"} for product rule,
#'     returns a full grid (default).
#'      \code{"sparse"} for combination technique,
#'       leads to a regular sparse grid.}
#'   \item{\code{level.trans}}{Logical variable denotes either to take the levels as number of grid points (FALSE = default) or to transform in that manner that number of grid points = 2^(levels-1) (TRUE). See, code{\link[mvQuad]{createNIGrid}}, for details.}
#' }
#'  \code{quadrature} is activated when \code{crt_method = "quadrature"} in
#'       any of the parent functions (for example, \code{\link{bayes}}).
#'
#' @examples
#' crt.bayes.control()
#' crt.bayes.control(cubature = list(tol = 1e-4))
#' crt.bayes.control(quadrature = list(level = 4))
#' @return A list  of two lists each contains the  control parameters for \code{\link[cubature]{hcubature}} and \code{\link[mvQuad]{createNIGrid}}, respectively.
#' @export
crt.bayes.control <- function(cubature = list(tol = 1e-5,
                                              maxEval = 50000,
                                              absError = 0),
                              quadrature = list(type = "GLe",
                                                level = 6,
                                                ndConstruction = "product",
                                                level.trans = FALSE
                                                )){
  cubature_out <- do.call(control.cubature, cubature)
  if (is.null(cubature$tol))
    cubature_out$tol <- 1e-5
  if (is.null(cubature$maxEval))
    cubature_out$maxEval <- 50000
  if (is.null(cubature$absError))
    cubature_out$absError <- 0

  ## qudrature part
  quadrature_out <- do.call(control.quadrature, quadrature)
  if (is.null(quadrature$type))
    quadrature_out$type <- "GLe"
  if (is.null(quadrature$level))
    quadrature_out$level <- 6
  if (is.null(quadrature$ndConstruction))
    quadrature_out$ndConstruction <- "product"
  if (is.null(quadrature$level.trans))
    quadrature_out$level.trans <- FALSE

  return(list(cubature = cubature_out, quadrature = quadrature_out))
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
#' @importFrom mvQuad createNIGrid
#' @importFrom mvQuad rescale
#' @importFrom mvQuad quadrature
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
  type <- arg$type
  #if (!arg$is.only.w)
  npred <- length(arg$lx)
  #else
  #   npred <- NA
  ## number of parameters
  #npar <- arg$npar
  if (!(type %in% c("D", "DPA", "DPM", "multiple", "D_LLTM")))
    stop("bug: 'type' must be  'D' or 'DPM' or 'DPM' or 'multiple' or 'D_LLTM' in 'iterate.ICAB")
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
  if (!arg$is.only.w){
    if (ICA.control$sym)
      x_id <- 1:floor(arg$k/2) else
        x_id <- 1:(arg$k * npred)
      if (!ICA.control$equal_weight)
        w_id <- (x_id[length(x_id)] + 1):length(arg$ld) else
          w_id <- NA
  }else{
    w_id <- 1:length(arg$ld)
    x_id <- NA
  }


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
                   Calculate_Cost = Calculate_Cost_bayes,
                   is.only.w = arg$is.only.w)

  vertices_outer <- make_vertices(lower = arg$lx, upper = arg$ux)
  sens_varlist <-list(npred = npred,
                      # plot_3d = "lattice",
                      npar = arg$npar,
                      fimfunc_sens = arg$FIM_sens,
                      Psi_x_bayes  = arg$Psi_funcs$Psi_x_bayes,
                      Psi_xy_bayes  = arg$Psi_funcs$Psi_xy_bayes,
                      crfunc = arg$crfunc,
                      vertices_outer = vertices_outer,
                      is.only.w = arg$is.only.w)
  ########################################################################################*
  #cat("iterate ", get(".Random.seed")[2], "\n")
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
    prev_time <- 0
    ############################################## Initialization for ICA
    InitialCountries <- GenerateNewCountry(NumOfCountries = ICA.control$ncount,
                                           lower = arg$ld,
                                           upper = arg$ud,
                                           sym = ICA.control$sym,
                                           w_id = w_id,
                                           x_id = x_id,
                                           npred= npred,
                                           equal_weight = ICA.control$equal_weight,
                                           is.only.w = arg$is.only.w)
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
    if (exists(".Random.seed")){
      GlobalSeed <- get(".Random.seed", envir = .GlobalEnv)
      #if you call directly from iterate and not bayes!
      on.exit(assign(".Random.seed", GlobalSeed, envir = .GlobalEnv))
    }
    msg <- object$best$msg
    prev_iter <- length(evol) ##previous number of iterationst
    maxiter <- iter + prev_iter
    totaliter <- prev_iter
    mean_cost <- sapply(1:(totaliter), FUN = function(j) evol[[j]]$mean_cost)
    min_cost <- sapply(1:(totaliter), FUN = function(j) evol[[j]]$min_cost)
    Empires <- object$empires
    prev_time <- arg$time

    check_counter <- arg$updating$check_counter
    total_nfeval <- object$alg$nfeval
    total_nlocal <-  object$alg$nlocal
    total_nrevol <-object$alg$nrevol
    total_nimprove <-object$alg$nimprove

    imp_cost <- round(sapply(object$empires, "[[", "ImperialistCost"), 12)
    best_imp_id<- which.min(imp_cost)
    revol_rate <-  arg$updating$revol_rate
    ##updating the random seed

    if (!is.null(ICA.control$rseed)){
      do.call("RNGkind",as.list(arg$updating$oldRNGkind))  ## must be first!
      assign(".Random.seed", arg$updating$oldseed , .GlobalEnv)
    }
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
        # if (ii == 4){
        #   cat(Empires[[ii]]$ImperialistPosition)
        #   return(NULL)
        # }

        #cat("\nbefore local: empire",ii, " des:",  Empires[[ii]]$ImperialistPosition)
        LocalSearch_res <- LocalSearch (TheEmpire =  Empires[[ii]],
                                        lower = arg$ld,
                                        upper = arg$ud,
                                        l = ICA.control$l,
                                        fixed_arg = fixed_arg)


        Empires[[ii]] <- LocalSearch_res$TheEmpire



        # cat("\nafter local: empire",ii, " des:",  Empires[[ii]]$ImperialistPosition)

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
      #cat("\nafter assim: empire",ii, " des:",  Empires[[ii]]$ImperialistPosition)
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

      #cat("\nafter revol: empire",ii, " des:",  Empires[[ii]]$ImperialistPosition)
      ############################################################*
      Empires[[ii]] <- PossesEmpire(TheEmpire = Empires[[ii]])
      # cat("\nafter possession: empire",ii, " des:",  Empires[[ii]]$ImperialistPosition)
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
    if (type %in% c("D", "DPA", "DPM", "multiple", "D_LLTM")){
      min_cost[totaliter] <-   min(imp_cost)
      mean_cost[totaliter] <-  mean(imp_cost)
    }else
      stop("Bug: check the type in iterate.bayes")
    best_imp_id <- which.min(imp_cost) ## which list contain the best imp
    if (!ICA.control$equal_weight)
      w <- Empires[[best_imp_id]]$ImperialistPosition[, w_id] else
        w <- w_equal
    if (!arg$is.only.w)
      x <- Empires[[best_imp_id]]$ImperialistPosition[, x_id] else
        x <- arg$only_w_varlist$x

    if (ICA.control$sym){
      x_w <- ICA_extract_x_w(x = x, w = w, sym_point = ICA.control$sym_point)
      x <- x_w$x
      w <- x_w$w
    }


    ##sort Point
    if (!arg$is.only.w){
      if (npred == 1){
        w <- w[order(x)]
        x <- sort(x)
      }
    }
    ############################################################################*

    ################################################################ print trace

    if (ICA.control$trace){
      if (!arg$is.only.w)
        cat("\nIteration:", totaliter, "\nDesign Points:\n", x, "\nWeights: \n", w,
            "\nCriterion value: ", min_cost[totaliter],
            "\nTotal number of function evaluations:", total_nfeval, "\nTotal number of successful local search moves:", total_nlocal,
            "\nTotal number of successful revolution moves:", total_nrevol, "\n") else
              cat("\nIteration:", totaliter, "\nWeights: \n", w,
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
      check_counter <- 0
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

  object$arg$time <- proc.time() - arg$time_start + prev_time
  ## arg has a list named update as well
  ###### end of saving
  ##############################################################################*
  return(object)

}


######################################################################################################*
######################################################################################################*



######################################################################################################*
######################################################################################################*
