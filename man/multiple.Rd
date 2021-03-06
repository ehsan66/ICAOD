% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3-UserMinimaxFunctions.R
\name{multiple}
\alias{multiple}
\title{Locally Multiple Objective Optimal Designs for the 4-Parameter Hill Model}
\usage{
multiple(
  minDose,
  maxDose,
  iter,
  k,
  inipars,
  Hill_par = TRUE,
  delta,
  lambda,
  fimfunc = NULL,
  ICA.control = list(),
  sens.control = list(),
  initial = NULL,
  tol = sqrt(.Machine$double.xmin),
  x = NULL
)
}
\arguments{
\item{minDose}{Minimum dose \eqn{D}. For the 4-parameter logistic model, i.e. when \code{Hill_par = FALSE}, it is the minimum of \eqn{log(D)}.}

\item{maxDose}{Maximum dose \eqn{D}. For the 4-parameter logistic model, i.e. when \code{Hill_par = FALSE}, it is the maximum of \eqn{log(D)}.}

\item{iter}{Maximum number of iterations.}

\item{k}{Number of design points. Must be at least equal to the number of model parameters to avoid singularity of the FIM.}

\item{inipars}{A vector of initial estimates for the vector of parameters  \eqn{(a, b, c, d)}.
For the 4-parameter logistic model, i.e. when \code{Hill_par = FALSE},
it is  a vector of initial estimates for \eqn{(\theta_1, \theta_2,\theta_3, \theta_4)}.}

\item{Hill_par}{Hill model parameterization? Defaults to \code{TRUE}.}

\item{delta}{Predetermined meaningful value of the minimum effective dose MED.
When \eqn{\delta < 0 }, then \eqn{\theta_2 > 0} or when \eqn{\delta > 0}, then \eqn{\theta_2 < 0}.}

\item{lambda}{A vector of relative importance of each of the three criteria,
i.e. \eqn{\lambda = (\lambda_1, \lambda_2, \lambda_3)}.
 Here \eqn{0 < \lambda_i < 1} and  s \eqn{\sum \lambda_i = 1}.}

\item{fimfunc}{A function. Returns the FIM as a \code{matrix}. Required when \code{formula} is missing. See 'Details' of \code{\link{minimax}}.}

\item{ICA.control}{ICA control parameters. For details, see \code{\link{ICA.control}}.}

\item{sens.control}{Control Parameters for Calculating the ELB. For details, see \code{\link{sens.control}}.}

\item{initial}{A matrix of the  initial design points and weights that will be inserted into the initial solutions (countries) of the algorithm.
Every row is a design, i.e.  a concatenation of \code{x} and \code{w}. Will be coerced to a \code{matrix} if necessary.  See 'Details' of \code{\link{minimax}}.}

\item{tol}{Tolerance for finding the general inverse of the Fisher information matrix. Defaults to \code{.Machine$double.xmin}.}

\item{x}{A vector of candidate design (support) points.
When is not set to \code{NULL} (default),
 the algorithm only finds the optimal weights for the candidate points in  \code{x}.
   Should be set when the user has a finite number of candidate design points  and the purpose
   is to find the optimal weight for each of them (when zero, they will be excluded from the design).
For design points with more than one dimension, see 'Details' of \code{\link{sensminimax}}.}
}
\value{
an object of class \code{minimax} that is a list including three sub-lists:
\describe{
  \item{\code{arg}}{A list of design and algorithm parameters.}
  \item{\code{evol}}{A list of length equal to the number of iterations that stores
   the information about the best design (design with least criterion value)
    of each iteration. \code{evol[[iter]]} contains:
    \tabular{lll}{
      \code{iter}                   \tab      \tab Iteration number.\cr
      \code{x}                      \tab      \tab Design points. \cr
      \code{w}                      \tab      \tab Design weights. \cr
      \code{min_cost}               \tab      \tab Value of the criterion for the best imperialist (design).  \cr
      \code{mean_cost}              \tab      \tab Mean of the criterion values of all the imperialists. \cr
      \code{sens}                   \tab      \tab An object of class \code{'sensminimax'}. See below. \cr
      \code{param}                  \tab      \tab Vector of parameters.\cr
    }
  }

  \item{\code{empires}}{A list of all the  empires of the last iteration.}
  \item{\code{alg}}{A list with following information:
    \tabular{lll}{
      \code{nfeval}           \tab      \tab Number of function evaluations.  It does not count the function evaluations from checking the general equivalence theorem.\cr
      \code{nlocal}           \tab      \tab Number of successful local searches. \cr
      \code{nrevol}           \tab      \tab Number of successful revolutions. \cr
      \code{nimprove}         \tab      \tab Number of successful movements toward the imperialists in the assimilation step. \cr
      \code{convergence}      \tab      \tab Stopped by \code{'maxiter'} or \code{'equivalence'}?\cr
    }
  }
  \item{\code{method}}{A type of optimal designs used.}
  \item{\code{design}}{Design points and weights at the final iteration.}
  \item{\code{out}}{A data frame of design points, weights, value of the criterion for the best imperialist (min_cost), and Mean of the criterion values of all the imperialistsat each iteration (mean_cost).}
}

The list \code{sens} contains information about the design verification by the general equivalence theorem. See \code{sensminimax} for more details.
It is given  every \code{ICA.control$checkfreq} iterations
and also the last iteration if   \code{ICA.control$checkfreq >= 0}. Otherwise, \code{NULL}.

 \code{param} is a vector of parameters that is the global minimum of
  the minimax criterion or the global maximum of the standardized maximin criterion over the parameter space, given  the current \code{x}, \code{w}.
}
\description{
The 4-parameter Hill model is of the form
 \deqn{f(D) = c + \frac{(d-c)(\frac{D}{a})^b}{1+(\frac{D}{a})^b} + \epsilon,}{
 f(D) = c + (d-c)(D/a)^b/(1 + (D/a)^b) + \epsilon,}
where \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)},
 \eqn{D} is the dose level and the predictor,
\eqn{a} is the ED50,
 \eqn{d} is the upper limit of response,
  \eqn{c} is the lower limit of response and
   \eqn{b} denotes the Hill constant that
 control the flexibility in the slope of the response curve.\cr
 Sometimes, the Hill model is re-parameterized and written as
 \deqn{f(x) = \frac{\theta_1}{1 + exp(\theta_2 x + \theta_3)} + \theta_4,}{
 f(x)= \theta_1/(1 + exp(\theta_2*x + \theta_3)) + \theta_4,}
  where \eqn{\theta_1 = d - c}, \eqn{\theta_2 = - b},
  \eqn{\theta_3 = b\log(a)}{\theta_3 = b*log(a)}, \eqn{\theta_4 = c}, \eqn{\theta_1 > 0},
  \eqn{\theta_2 \neq 0}{\theta_2 not equal to 0}, and \eqn{-\infty < ED50 < \infty},
  where \eqn{x = log(D) \in [-M, M]}{x = log(D) belongs to [-M, M]}
  for some sufficiently large value of \eqn{M}.
  The new form is sometimes called  4-parameter logistic model.\cr
 The function \code{multiple} finds locally multiple-objective optimal designs for estimating the model parameters, the ED50, and the MED, simultaneously.
   For more details, see Hyun and  Wong (2015).
}
\details{
When \eqn{\lambda_1 > 0}, then the number of support points \code{k}
  must at least be four to avoid singularity of the Fisher information matrix.

One can adjust the tuning parameters in \code{\link{ICA.control}} to set a stopping rule
based on the general equivalence theorem. See 'Examples' below.
}
\note{
This function is NOT appropriate for  finding  c-optimal designs for estimating 'MED' or 'ED50' (single objective optimal designs)
 and  the results may not be stable.
 The reason is that for the c-optimal criterion
 the generalized inverse of the Fisher information matrix
  is not stable and depends
 on the tolerance value (\code{tol}).
}
\examples{
# All the examples are available in Hyun and Wong (2015)

#################################
#  4-parameter logistic model
# Example 1, Table 3
#################################
lam <- c(0.05, 0.05, .90)
# Initial estimates are derived from Table 1
# See how the stopping rules are set via 'stop_rul', checkfreq' and 'stoptol'
Theta1 <- c(1.563, 1.790, 8.442, 0.137)
res1 <- multiple(minDose = log(.001), maxDose = log(1000),
                 inipars = Theta1, k = 4, lambda = lam, delta = -1,
                 Hill_par = FALSE,
                 iter = 1,
                 ICA.control = list(rseed = 1366, ncount = 100,
                                    stop_rule = "equivalence",
                                    checkfreq = 100, stoptol = .95))
\dontrun{
res1 <- update(res1, 1000)
# stops at iteration 101
}

#################################
#  4-parameter Hill model
#################################
## initial estimates for the parameters of Hill model:
a <- 0.008949  # ED50
b <- -1.79 # Hill constant
c <- 0.137 # lower limit
d <- 1.7 # upper limit
# D belongs to c(.001, 1000) ## dose in mg
## the vector of Hill parameters are now c(a, b, c, d)
\dontrun{
res2 <- multiple(minDose = .001, maxDose = 1000,
                 inipars =  c(a, b, c, d),
                 Hill_par = TRUE, k = 4, lambda = lam,
                 delta = -1, iter = 1000,
                 ICA.control = list(rseed = 1366, ncount = 100,
                                    stop_rule = "equivalence",
                                    checkfreq = 100, stoptol = .95))
# stops at iteration 100
}



# use x argument to provide fix number of  dose levels.
# In this case, the optimization is only over weights
\dontrun{
res3 <- multiple(minDose = log(.001), maxDose = log(1000),
                 inipars = Theta1, k = 4, lambda = lam, delta = -1,
                 iter = 300,
                 Hill_par = FALSE,
                 x = c(-6.90, -4.66, -3.93, 3.61),
                 ICA.control = list(rseed = 1366))
res3$evol[[300]]$w
# if the user provide the desugn points via x, there is no guarantee
#   that the resulted design is optimal. It only provides the optimal weights given
#   the x points of the design.
plot(res3)
}

}
\references{
Hyun, S. W., and Wong, W. K. (2015). Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. The international journal of biostatistics, 11(2), 253-271.
}
\seealso{
\code{\link{sensmultiple}}
}
