% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8-UserPriors.R
\name{uniform}
\alias{uniform}
\title{Assume A Multivariate Uniform Prior Distribution for The Model Parameters}
\usage{
uniform(lower, upper)
}
\arguments{
\item{lower}{A vector of lower bounds  for the model parameters.}

\item{upper}{A vector of upper bounds  for the model parameters.}
}
\value{
An object of class \code{cprior} that is a  list with the following components:
\itemize{
 \item{\code{fn}: }{prior distribution as an R \code{function} with argument \code{param} that is the vector of the unknown parameters. See below.}
 \item{\code{npar}: }{Number of unknown parameters and is equal to the length of \code{param}}.
 \item{\code{lower}: }{Argument \code{lower}. It has the same length as \code{param}}.
 \item{\code{upper}: }{Argument \code{lower}. It has the same length as \code{param}}.
}
The list will be passed to the argument \code{prior} of the function \code{\link{bayes}}.
 The order of the argument \code{param} in \code{fn} has the same order as the argument \code{parvars} when the model is specified by a formula.
Otherwise, it is equal to the argument \code{param} in the function \code{fimfunc}.
}
\description{
Creates independent uniform prior distributions for the unknown model parameters as an object of class \code{cprior}.
}
\note{
The order of the argument \code{param} in \code{fn} has the same order as the argument \code{parvars} when the model is specified by a formula.
Otherwise, it is the same as the argument \code{param} in the function \code{fimfunc}.
}
\examples{
uniform(lower =  c(-3, .1), upper = c(3, 2))
}
\seealso{
\code{\link{bayes}} \code{\link{sensbayes}}
}
