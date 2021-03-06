% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{FIM_logistic_2pred}
\alias{FIM_logistic_2pred}
\title{Fisher Information Matrix for the Logistic Model with Two Predictors}
\usage{
FIM_logistic_2pred(x1, x2, w, param)
}
\arguments{
\item{x1}{Vector of design points (for first predictor).}

\item{x2}{Vector of design points (for second predictor).}

\item{w}{Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.}

\item{param}{Vector of values for the model parameters \code{c(b0, b1, b2, b3)}.}
}
\value{
Fisher information matrix.
}
\description{
It provides the cpp function for FIM for the following model:\cr
  \code{~exp(b0+ b1 * x1 + b2 * x2 + b3 * x1 * x2)/(1 + exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))}.
}
