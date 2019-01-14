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


res <- vector("list", 3)
############################################################################*
# kinetics of the catalytic dehydra- tion of n-hexyl alcohol -----
############################################################################*
# Duarte et al. (2016), Table 8

#################################
# uniform
#################################
res[[1]] <- bayes(formula =  ~b3*x1/(1+b1*x1 + b2*x2),
                  predvars = c("x1", "x2"),
                  parvars = c("b1", "b2", "b3"),
                  lx = rep(0, 2), ux = rep(2, 2),
                  prior = uniform(lower = c(1.9, 9.2, 1.14),
                                  upper = c(3.9, 15.2, 2.34)),
                  k = 3,
                  ICA.control = list(rseed = 1366),
                  iter = 150)
plot(res[[1]])
res[[1]]$arg$time
#################################
# normal
#################################
res[[2]] <- bayes(formula =  ~b3*x1/(1+b1*x1 + b2*x2),
                  predvars = c("x1", "x2"),
                  parvars = c("b1", "b2", "b3"),
                  lx = rep(0, 2), ux = rep(2, 2),
                  prior = normal(mu = c(2.9, 12.2, 1.74),
                                 sigma = diag(c((1/3)^2, 1^2, 0.2^2)),
                                 lower = c(1.9, 9.2, 1.14),
                                 upper = c(3.9, 15.2, 2.34)),
                  k = 3,
                  ICA.control = list(rseed = 1366),
                  iter = 150)
res[[2]]$arg$time

# less value for maxEval
res[[3]] <- bayes(formula =  ~b3*x1/(1+b1*x1 + b2*x2),
                  predvars = c("x1", "x2"),
                  parvars = c("b1", "b2", "b3"),
                  lx = rep(0, 2), ux = rep(2, 2),
                  prior = normal(mu = c(2.9, 12.2, 1.74),
                                 sigma = diag(c((1/3)^2, 1^2, 0.2^2)),
                                 lower = c(1.9, 9.2, 1.14),
                                 upper = c(3.9, 15.2, 2.34)),
                  k = 3,
                  crt.bayes.control = list(cubature = list(maxEval=2000, tol = 1e-5)),
                  ICA.control = list(rseed = 1366),
                  iter = 150)
res[[3]]$arg$time
