######################################################################################################*
######################################################################################################*
#'@title Calculates Relative Efficiency for Locally Optimal Designs
#' @description
#' Given a vector of initial estimates for the parameters, this function calculates the D-and PA- efficiency of a design \eqn{\xi_1} with respect to a design \eqn{\xi_2}.
#' Usually, \eqn{\xi_2} is an  optimal design.
#'
#'
#' @details For a known \eqn{\theta_0}, relative D-efficiency is
#' \deqn{exp(\frac{log|M(\xi, \theta_0)| - log|M(\xi^*, \theta_0)|}{npar})}{
#' exp\{(log|M(\xi, \theta_0)| - log|M(\xi*, \theta_0)|)/npar\}.}
#' The relative P-efficiency is
#' \deqn{\exp(\log(\sum_{i=1}^k w_ip(x_i, \theta_0) - \log(\sum_{i=1}^k w^*_ip(x^*_i, \theta_0))}{
#' exp(log (\sum w_i p(x_i, \theta_0) - log(\sum w*_i p(x*_i, \theta_0)),
#' }
#' where \eqn{x^*}{x*} and \eqn{w^*}{w*} are the support points and the corresponding weights of the optimal design, respectively.
#'
#' The argument  \code{x} is the vector of design points.
#'  For design points with more than one dimension (the models with more than one predictors),
#'    it is a concatenation of the design points, but \strong{dimension-wise}.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  a two-point optimal design has the following points:
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     Then, the argument \code{x} is equal to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'
#' @export
#' @inheritParams locallycomp
#' @param x Vector of design (support) points of \eqn{\xi_1}. See 'Details' of \code{\link{leff}}.
#' @param w Vector of corresponding design weights for \code{x}.
#' @param xopt Vector of design (support) points of the optimal design (\eqn{\xi_2}). Similar to \code{x}.
#' @param wopt Vector of corresponding design weights for \code{xopt}.
#' @param type A character. \code{"D"} denotes the D-efficiency and \code{"PA"} denotes the average P-efficiency.
#' @return A value between 0 and 1.
#' @references McGree, J. M., Eccleston, J. A., and Duffull, S. B. (2008). Compound optimal design criteria for nonlinear models. Journal of Biopharmaceutical Statistics, 18(4), 646-661.
#' @example inst/examples/leff_examples.R
leff <- function(formula,
                 predvars,
                 parvars,
                 family = gaussian(),
                 inipars,
                 type = c("D", "PA"),
                 fimfunc = NULL,
                 xopt, wopt, x, w,
                 npar = length(inipars),
                 prob = NULL){
  ## bayesian D-efficiency
  ### relative efficieny of x with respect to xopt
  ## if the user use small p

  if (!(type[1] %in% c("D", "PA")))
    stop("'type' must be 'D' or 'PA'")
  npred <- length(x)/length(w)

  fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =rep(0, npred), ux = rep(1, npred), iter = 1, k = 1,
                                       paramvectorized = FALSE,
                                       prior = NULL, x = NULL)


  if (!missing(formula)){
    if (length(inipars) != length(parvars))
      stop("length of 'inipars' is not equal to the length of 'parvars'")
  }

  if(missing(formula)){
    # to handle ...
    fimfunc2 <- function(x, w, param)
      fimfunc(x = x, w = w, param = param)
    #fimfunc(x = x, w = w, param = param,...)
  } else{
    fimfunc2 <- fimfunc_formula$fimfunc_formula ## can be vectorized with respect to parameters!
  }
  if (type[1] == "D")
    releff <- (det2(fimfunc2(x = x, w = w, param = inipars))/
                 det2(fimfunc2(x = xopt, w = wopt, param = inipars)))^(1/npar)

  if (type[1] == "PA"){

    if (is.null(prob))
      stop("'prob' must be given for P-optimality")

    if (is.formula(prob)){
      prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
    }else{
      if (!is.function(prob))
        stop("'prob' must be either a function or a formula")
      if (!formalArgs(prob) %in% c("x", "param"))
        stop("arguments of 'prob' must be 'x' and 'param'")
    }
    releff <- sum(w * prob(x = x, param = inipars))/sum(wopt * prob(x = xopt, param = inipars))
  }


  return(releff)
}

######################################################################################################*
######################################################################################################*
#' @title  Calculates Relative Efficiency for Bayesian Optimal Designs
#' @description
#' Given a prior distribution for the parameters, this function calculates the Bayesian D-and PA- efficiency of a design \eqn{\xi_1} with respect to a design \eqn{\xi_2}.
#' Usually, \eqn{\xi_2} is an optimal design.
#' This function is especially useful for investigating the robustness of Bayesian optimal designs under different prior distributions (See 'Examples').
#'
#' @inheritParams leff
#' @inheritParams bayes
#' @details
#' See Masoudi et al. (2018) for formula details (the paper is under review and will be updated as soon as accepted).
#'
#' The argument  \code{x} is the vector of design points.
#'  For design points with more than one dimension (the models with more than one predictors),
#'    it is a concatenation of the design points, but \strong{dimension-wise}.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  a two-point optimal design has the following points:
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     Then, the argument \code{x} is equal to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'
#'
#' @export
#' @example inst/examples/beff_examples.R
beff <- function(formula,
                 predvars,
                 parvars,
                 family = gaussian(),
                 prior,
                 fimfunc = NULL,
                 xopt, wopt, x, w,
                 crt.bayes.control = list(),
                 npar = NULL,
                 type = c("D", "PA"),
                 prob = NULL){
  ## bayesian D-efficiency
  ### relative efficieny of x with respect to xopt
  if (!(type[1] %in% c("D", "PA")))
    stop("'type' must be 'D' or 'PA'")
  if (is.null(npar)){
    if (!missing(formula))
      npar <- length(parvars)
    else
      npar <- prior$npar
  }

  ## only to pass the check_common_eargs
  npred <- length(x)/length(w)
  fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =rep(0, npred), ux = rep(1,npred),
                                       iter = 1, k = length(w),
                                       paramvectorized = FALSE, prior = prior,
                                       x = NULL)
  if(missing(formula)){
    # to handle ...
    fimfunc2 <- function(x, w, param){
      fimfunc(x = x, w = w, param = param)
      #fimfunc(x = x, w = w, param = param,...)
    }
  } else{

    #if (length(prior$lower) != length(parvars))
    if (length(prior$lower) != fimfunc_formula$num_unknown_param)
      stop("length of 'prior$lower' is not equal to the number of unknown (not fixed) parameters")
    # fim_localdes <- fimfunc_formula$fimfunc_formula
    fimfunc2 <- fimfunc_formula$fimfunc_formula ## can be vectorized with respect to parameters!
  }
  crt.bayes.control <- do.call("crt.bayes.control", crt.bayes.control)

  if (type[1] == "D")
    cr_integrand <- function(param, x, w){
      # bcrfunc1 <- apply(param, 2,
      #                   FUN = function(col_par)-det2(fimfunc2(x = x, w = w, param = col_par), logarithm = TRUE)) * prior$fn(t(param))
      bcrfunc1 <- apply(param, 2,
                        FUN = function(col_par)-det2(fimfunc2(x = x, w = w, param = col_par), logarithm = TRUE)) * prior$fn(t(param))

      dim(bcrfunc1) <- c(1, length(bcrfunc1))
      return(bcrfunc1)
    }
  if (type[1] == "PA"){
    if (is.formula(prob)){
      prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
    }else{
      if (!is.function(prob))
        stop("'prob' must be either a function or a formula when type is 'PA'")
      if (!formalArgs(prob) %in% c("x", "param"))
        stop("arguments of 'prob' must be 'x' and 'param'")
    }
    cr_integrand <- function(param, x, w){
      bcrfunc1 <- (apply(param, 2, function(col_par) -log(sum(w * prob(x = x, param = col_par))))) * prior$fn(t(param))
      dim(bcrfunc1) <- c(1, length(bcrfunc1))
      return(bcrfunc1)
    }

  }
  crfunc_bayesian  <- function(x, w, maxEval, tol) {

    out <- hcubature(f = cr_integrand, lowerLimit = prior$lower,
                     upperLimit = prior$upper,
                     vectorInterface = TRUE,
                     x = x, w = w, tol = tol, maxEval = maxEval)

    val <- out$integral
    return(list(val = val, fneval = out$functionEvaluations))
  }

  if (type[1] == "D")
    releff<- exp((crfunc_bayesian(x = xopt, w = wopt, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val-
                    crfunc_bayesian(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)/npar)
  # if (type[1] == "PA")
  #   releff<-  (crfunc_bayesian(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val/
  #                   crfunc_bayesian(x = xopt, w = wopt, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)
  if (type[1] == "PA")
    releff <- exp(crfunc_bayesian(x = xopt, w = wopt, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val - crfunc_bayesian(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)

  # releff<- crfunc_bayesian_D(x = xopt, w = wopt, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val/
  #                crfunc_bayesian_D(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val

  return(releff)
}
######################################################################################################*
######################################################################################################*

######################################################################################################*
######################################################################################################*

######################################################################################################*
######################################################################################################*
