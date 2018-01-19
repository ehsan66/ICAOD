######################################################################################################*
######################################################################################################*
#' Multivariate Normal Prior Distribution for Model Parameters
#'
#' Creates a multivariate normal prior distribution for the unknown parameters as an object of class \code{cprior}.
#'
#' @param mu  A vector of length of parameters, representing the mean value.
#' @param sigma A symmetric positive-definite matrix representing the variance-covariance matrix of the distribution.
#' @param lower A vector of lower bounds  for the unknown parameters.
#' @param upper A vector of upper bounds  for the unknown parameters.
#' @return
#' An object of class \code{cprior} that is a  list with components:
#' \itemize{
#'  \item{fn: }{prior distribution as an R \code{function} with argument \code{param} that is the vector of the unknown parameters. See below.}
#'  \item{npar: }{Number of unknown parameters and is equal to the length of \code{param}}.
#'  \item{lower: }{Argument \code{lower}. It has the same length as \code{param}}.
#'  \item{upper: }{Argument \code{lower}. It has the same length as \code{param}}.
#' }
#' The list will be passed to the argument \code{prior} of the function \code{\link{bayes}}.
#'  The order of the argument \code{param} in \code{fn} has the same order as the argument \code{parvars} when the model is specified by a formula.
#' Otherwise, it is the same as the argument \code{param} in the function \code{fimfunc}.
#' @export
#' @seealso \code{\link{bayes}} \code{\link{sensbayes}}
## @importFrom mnormt dmnorm
#' @examples normal(mu =  c(0, 1), sigma = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
#'   lower =  c(-3, .1), upper = c(3, 2))
normal <- function(mu, sigma, lower, upper){

  if (!is.matrix(sigma)){
    if (length(mu) != length(sigma))
      stop("check the input of 'normal' prior function. length of 'mu' is not equal to length of 'sigma'")
  }else
    if (ncol(sigma) != length(mu))
      stop("length of mean is not equal to the number of columns of covariance matrix")
  npar <- length(mu)
  prior_func  <- NA ## to define the variable in the global environment and avoid R CMD check Note
  prior_char <- paste("prior_func <- function(param){ \n out <- mnormt::dmnorm(x = param, mean = c(", paste(mu, collapse = ", "),
                      "), varcov = matrix(c(", paste(sigma, collapse = ", "), "), nrow =", npar , "))\n return(matrix(out, ncol = dim(param)[1]))}", sep = "")
  eval(parse(text = prior_char))
  return(structure(list(fn = prior_func, npar = npar, lower = lower, upper = upper), class =  "cprior"))
}





######################################################################################################*
######################################################################################################*
#' Multivariate Uniform Prior Distribution for Model Parameters
#'
#' Creates independent uniform prior distributions for the unknown model parameters as an object of class \code{cprior}.
#'
#' @inheritParams normal
#' @note The order of the argument \code{param} in \code{fn} has the same order as the argument \code{parvars} when the model is specified by a formula.
#' Otherwise, it is the same as the argument \code{param} in the function \code{fimfunc}.
#' @export
#' @inherit normal return
#' @seealso \code{\link{bayes}} \code{\link{sensbayes}}
#' @examples uniform(lower =  c(-3, .1), upper = c(3, 2))
uniform <- function(lower, upper){
  if (length(lower) != length(upper))
    stop(" length of 'lower' is not equal to the length of 'upper'")
  npar <- length(lower)
  prior_func  <- NA ## to define the variable in the global environment and avoid R CMD check Note
  prior_char <- paste("prior_func <- function(param){ return(matrix(c(",  1/prod((upper - lower)), "), nrow = dim(param)[1]))}", sep = "")
  eval(parse(text = prior_char))
  return(structure(list(fn = prior_func, npar = npar, lower = lower, upper = upper), class =  "cprior"))
}
######################################################################################################*
######################################################################################################*
#' Multivariate Skewed Normal Prior Distribution for Model Parameters
#'
#' Creates a multivariate skewed normal prior distribution for the unknown parameters as an object of class \code{cprior}.
#'
#' @inheritParams normal
#' @param xi A numeric vector of length \code{d=length(alpha)} representing the location parameter of the distribution. See 'Background' in \code{\link[sn]{dmsn}}.
#' @param Omega A symmetric positive-definite matrix of dimension \code{(d,d)}. See 'Background' in \code{\link[sn]{dmsn}}.
#' @param alpha A numeric vector which regulates the slant of the density. See 'Background' in \code{\link[sn]{dmsn}}.
#' @export
#' @inherit normal return
#' @importFrom sn dmsn
#' @seealso \code{\link{bayes}} \code{\link{sensbayes}}
#' @examples skewnormal(xi = c(0, 1),
#'  Omega = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
#'   alpha = c(1, 0), lower =  c(-3, .1), upper = c(3, 2))
skewnormal <- function(xi, Omega, alpha, lower, upper){
  npar <- length(alpha)
  prior_func  <- NA ## to define the variable in the global environment and avoid R CMD check Note
  prior_char <- paste("prior_func <- function(param){\n  out <- sn::dmsn(x = param, xi = c(", paste(xi, collapse = ", "),
                      "), Omega = matrix(c(", paste(Omega, collapse = ", "), "), nrow = ", length(xi),
                      "), alpha = c(", paste(alpha, collapse = ", "),
                      "))\n return(matrix(out, ncol = dim(param)[1]))\n}", sep = "")
  eval(parse(text = prior_char))
  return(structure(list(fn = prior_func, npar = npar, lower = lower, upper = upper), class =  "cprior"))
}
######################################################################################################*
######################################################################################################*
#' Multivariate Student's t Prior Distribution for Model Parameters
#'
#' Creates the prior distribution for the parameters as an object of class \code{cprior}.
#'
#' @inheritParams normal
#' @param mean  A vector of length \code{d=ncol(S)}, representing the location parameter (equal to the mean vector when \code{df>1}). See 'Arguments' in \code{\link[mnormt]{dmt}}.
#' @param S A symmetric positive-definite matrix representing the scale matrix of the distribution, such that \code{S*df/(df-2)} is the variance-covariance matrix when \code{df>2}. See 'Arguments' in \code{\link[mnormt]{dmt}}.
#' @param df Degrees of freedom; it must be a positive integer. See 'Arguments' in \code{\link[mnormt]{dmt}}.
#' @export
#' @importFrom mnormt dmt
#' @inherit normal return
#' @seealso \code{\link{bayes}} \code{\link{sensbayes}}
#' @examples skewnormal(xi = c(0, 1),
#'  Omega = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
#'   alpha = c(1, 0), lower =  c(-3, .1), upper = c(3, 2))
student <- function(mean, S, df, lower, upper){
  npar <- ncol(S)
  prior_func  <- NA ## to define the variable in the global environment and avoid R CMD check Note
  prior_char <- paste("prior_func <- function(param){\n  out <- mnormt::dmt(x = param, mean = c(", paste(mean, collapse = ", "),
                      "), S = matrix(c(", paste(S, collapse = ", "), "), nrow = ", length(mean),
                      "), df = ", df,  ")\n return(matrix(out, ncol = dim(param)[1]))\n}", sep = "")
  eval(parse(text = prior_char))
  return(structure(list(fn = prior_func, npar = npar, lower = lower, upper = upper), class =  "cprior"))
}
######################################################################################################*
######################################################################################################*

######################################################################################################*
######################################################################################################*

