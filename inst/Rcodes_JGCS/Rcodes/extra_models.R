# to install and load the package ICAOD ----
install_and_load <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

install_and_load("ICAOD")

# This code produces Bayesian D-optimal designs for some extra models

#################################*
# Weibull model: Uniform prior ----
################################*
# see Dette, H., & Pepelyshev, A. (2008).
# Efficient experimental designs for sigmoidal growth models.
# Journal of statistical planning and inference, 138(1), 2-17.

## See how we fixed some parameters (a, b and h)
res19 <- bayes(formula = ~a - b * exp(-lambda * t ^h),
               predvars = c("t"),
               parvars = c("a=1", "b=1", "lambda", "h=1"),
               lx = .00001, ux = 20,
               prior = uniform(.5, 2.5), k = 5, iter = 400,
               ICA.control = list(rseed = 1366))
plot(res19)
res19$arg$time

#################################*
# Weibull model: Normal prior ----
################################*
norm3 <- normal(mu = 1, sigma = .1, lower = .5, upper = 2.5)
res20 <- bayes(formula = ~a - b * exp(-lambda * t ^h),
               predvars = c("t"),
               parvars = c("a=1", "b=1", "lambda", "h=1"),
               lx = .00001, ux = 20, prior = norm3, k = 4, iter = 400,
               ICA.control = list(rseed = 1366))


plot(res20)

res20$arg$time
#################################*
# Richards model: Normal prior ----
#################################*
norm4 <- normal(mu = c(1, 1), sigma = matrix(c(.2, 0.1, 0.1, .4), 2, 2),
                lower = c(.4, .4), upper = c(1.6, 1.6))

res21 <- bayes(formula = ~a/(1 + b * exp(-lambda*t))^h,
               predvars = c("t"),
               parvars = c("a=1", "b", "lambda", "h=1"),
               lx = .00001, ux = 10,
               prior = norm4,
               k = 5, iter = 400,
               ICA.control = list(rseed = 1366))

res21$arg$time
plot(res21) ## 88 seconds
#################################*
# Exponential model: Uniform prior ----
#################################*
res22 <- bayes(formula = ~a + exp(-b*x), predvars = "x",
               parvars = c("a = 1", "b"),
               lx = 0.0001, ux = 1,
               prior = uniform(lower = 1, upper = 20),
               iter = 400, k = 3,
               ICA.control= list(rseed = 1366))

plot(res22)
res22$arg$time

#################################*
# Power logistic model ----
#################################*
# See, Duarte, B. P., & Wong, W. K. (2014).
# A Semidefinite Programming based approach for finding
# Bayesian optimal designs for nonlinear models

res23 <- bayes(formula = ~1/(1 + exp(-b *(x - a)))^s, predvars = "x",
               parvars = c("a", "b", "s"),
               lx = -1, ux = 1,
               prior = uniform(lower = c(-.3, 6, .5), upper = c(.3, 8, 1)),
               k = 5, iter = 400)

plot(res23)

res23$arg$time

############################################################################*
# A two-variable generalized linear model with a gamma distributed response ----
############################################################################*
lb <- c(.5, 0, 0, 0, 0, 0)
ub <- c(2, 1, 1, 1, 1, 1)
formula1 <- ~beta0+beta1*x1+beta2*x2+beta3*x1^2+beta4*x2^2+beta5*x1*x2
res24 <- bayes(formula = formula1,
               predvars = c("x1", "x2"), parvars = paste("beta", 0:5, sep = ""),
               family = Gamma(),
               lx = rep(0, 2), ux = rep(1, 2),
               prior = uniform(lower = lb, upper = ub),
               k = 7,iter = 500, ICA.control = list(rseed = 1366))


res24$arg$time
plot(res24)
# 57.33 seconds
