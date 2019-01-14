# to install and load the package ICAOD ----
install_and_load <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

tun1 <- list(cubature = list(tol = 1e-5, maxEval = 10000))

install_and_load("ICAOD")



############################################################*
# Section S2.1: 4-parameter sigmoid Emax model ----
#############################################################*
# Theta1
res1 <- bayes(fimfunc = ICAOD::FIM_sig_emax,
              lx = .001, ux = 500, k = 4, iter = 200,
              prior = uniform(c(4, 11, 100, 5), c(5, 12, 105, 6)),
              ICA.control = list(rseed = 1366),
              crt_method = "quadrature",
              crt.bayes.control = list(quadrature = list(type = "GLe",
                                                         level = 2)))

plot(res1)


# Theta2
res2 <-  bayes(fimfunc = ICAOD::FIM_sig_emax,
               lx = .001, ux = 500, k = 4, iter = 200,
               prior = uniform(c(4, 11, 100, 5), c(6, 13, 115, 7)),
               ICA.control = list(rseed = 1366),
               crt_method = "quadrature",
               crt.bayes.control = list(quadrature = list(type = "GLe",
                                                          level = 2)))
plot(res2)



# Theta3
res3 <-  bayes(fimfunc = ICAOD::FIM_sig_emax,
               lx = .001, ux = 500, k = 5, iter = 400,
               prior = uniform(c(4, 11, 100, 5), c(8, 15, 130, 9)),
               ICA.control = list(rseed = 1366),
               crt_method = "quadrature",
               crt.bayes.control = list(quadrature = list(type = "GLe",
                                                          level = 4)))

plot(res3)

# Theta4
res4 <- bayes(fimfunc = ICAOD::FIM_sig_emax,
              lx = .001, ux = 500, k = 7, iter = 500,
              prior = uniform(c(4, 11, 100, 5), c(10, 18, 180, 11)),
              ICA.control = list(rseed = 13, ncount = 300, nimp = 30),
              crt_method = "quadrature",
              crt.bayes.control = list(quadrature = list(type = "GLe",
                                                         level = 10)))
############################################################################*
# Section S2.2: A two-variable generalized linear model with a gamma distributed response ----
############################################################################*
formula1 <- ~beta0+beta1*x1+beta2*x2+beta3*x1^2+beta4*x2^2+beta5*x1*x2


sigma <- diag(c(.25^2, rep(.16^2, 5)))
mu <- c(1.25, .5, .5, .5, .5, .5)
lb <- mu - 3 * sqrt(diag(sigma))
ub <- mu + 3 * sqrt(diag(sigma))
norm1 <- normal(sigma = sigma, mu = mu, lower = lb, upper = ub)

uni1 <- uniform(lower = c(.5, 0, 0, 0, 0, 0), upper = c(2, 1, 1, 1, 1, 1))
res_gam <- vector("list", 4)

res_gam[[1]]$crt_method <- "cubature"
res_gam[[1]]$crt_bayes_control = list(cubature = list(maxEval = 50000))
res_gam[[1]]$prior <- uni1


res_gam[[2]]$crt_method <- "quadrature"
res_gam[[2]]$crt_bayes_control = list(quadrature = list(level = 3))
res_gam[[2]]$prior <- uni1

res_gam[[3]]$crt_method <- "cubature"
res_gam[[3]]$crt_bayes_control = list(cubature = list(maxEval = 2000))
res_gam[[3]]$prior <- norm1

res_gam[[4]]$crt_method <- "quadrature"
res_gam[[4]]$crt_bayes_control = list(quadrature = list(level = 3))
res_gam[[4]]$prior <- norm1

for (i in 1:4)
  res_gam[[i]]$d <- bayes(formula = formula1,
                          predvars = c("x1", "x2"), parvars = paste("beta", 0:5, sep = ""),
                          family = Gamma(),
                          lx = rep(0, 2), ux = rep(1, 2),
                          prior = res_gam[[i]]$prior,
                          k = 7,iter = 500,
                          ICA.control = list(rseed = 1366, ncount = 40),
                          crt_method = res_gam[[i]]$crt_method,
                          crt.bayes.control = res_gam[[i]]$crt_bayes_control)


res_gam[[1]]$d$arg$time
res_gam[[2]]$d$arg$time
res_gam[[3]]$d$arg$time
res_gam[[4]]$d$arg$time
