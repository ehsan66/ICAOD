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


# This code produces all the optimal designs reported in the JCGS paper


#The Package ICAOD provides two interfaces for Fisher information matrices:
#1) fimfunc argument: user should obtain and write the FIM himself.
#2) using formula: based on the formula of response and the family, ICAOD creates the FIM automatically.

# The results from these two methods may slightly differ due to round off error.
# For the examples in the paper, we wrote each FIM in c++ and called the functions by Rcpp package.
# At the end of this script, there are some new models: Richards, Weibull, exponential, power logistic model and Aa two-variable generalized linear model with a gamma distributed response.


#############################################*
# Section 4.1: two-parameter logistic model ----
#############################################*
# define the prior distributon
uni <- uniform(lower =  c(-3, .1), upper = c(3, 2))
# xi_uni
res1 <- bayes(fimfunc = FIM_logistic,
              lx = -3, ux = 3,
              k =  5, iter = 1000, prior = uni,
              ICA.control = list(rseed = 1366))

# or let ICAOD create the FIM matrix: the results slightly differ due to round off error.
res1 <- bayes(formula = ~1/(1 + exp(-b *(x - a))),
              predvars = "x", parvars = c("a", "b"),
              family = binomial(),
              lx = -3, ux = 3,
              k =  5, iter = 1000, prior = uni,
              ICA.control = list(rseed = 1366))

#######*
# define the normal prior distributon
norm1 <- normal(mu =  c(0, 1),
                sigma = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                lower =  c(-3, .1), upper = c(3, 2))
#xi_norm1
res2 <- bayes(fimfunc = FIM_logistic, lx = -3, ux = 3, k =  4,
              iter = 1000, prior = norm1,
              ICA.control = list(rseed = 1366))
#######*
# defining the normal prior
norm2 <- normal(mu =  c(0, 1),
                sigma = matrix(c(1, 0, 0, .5), nrow = 2),
                lower =  c(-3, .1), upper = c(3, 2))
# xi_norm2
res3 <- bayes(fimfunc = FIM_logistic,
              lx = -3, ux = 3, k =  4, iter = 1000, prior = norm2,
              ICA.control = list(rseed = 1366))
#######*
# define the skewnorml prior
skew1 <- skewnormal(xi = c(0, 1),
                    Omega = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                    alpha = c(1, 0), lower =  c(-3, .1), upper = c(3, 2))
# xi_skew1
res4 <- bayes(fimfunc = FIM_logistic,
              lx = -3, ux = 3, k =  4, iter = 1000, prior = skew1,
              ICA.control = list(rseed = 1366, ncount = 80))


#######*
# define the prior
skew2 <- skewnormal(xi = c(0, 1),
                    Omega = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                    alpha = c(-1, 0), lower =  c(-3, .1), upper = c(3, 2))
# xi_skew2
res5 <- bayes(fimfunc = FIM_logistic,
              lx = -3, ux = 3, k =  4, iter = 1000, prior = skew2,
              ICA.control = list(rseed = 1366, ncount = 60))



# define the prior
stud <- student(mean =  c(0, 1), S   = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                df = 3, lower =  c(-3, .1), upper = c(3, 2))
# xi_t
res6 <- bayes(fimfunc = FIM_logistic,
              lx = -3, ux = 3, k =  5, iter = 1000, prior = stud,
              ICA.control = list(ncount = 200, rseed = 1366))



############################################################*
# Section 4.2: 4-parameter sigmoid Emax model ----
#############################################################*
# xi_Theta1
res7 <- bayes(fimfunc = ICAOD::FIM_sig_emax,
              lx = .001, ux = 500, k = 4, iter = 200,
              prior = uniform(c(4, 11, 100, 5), c(5, 12, 105, 6)),
              ICA.control = list(rseed = 1366))

# # or
# res7 <- bayes(formula = ~ theta1 + (theta2 - theta1)*(x^theta4)/(x^theta4 + theta3^theta4),
#               lx = .001, ux = 500, k = 4, iter = 200,
#               prior = uniform(c(4, 11, 100, 5), c(5, 12, 105, 6)),
#               ICA.control = list(rseed = 13))

# xi_Theta2
res8 <- bayes(fimfunc = ICAOD::FIM_sig_emax,
              lx = .001, ux = 500, k = 4, iter = 200,
              prior = uniform(c(4, 11, 100, 5), c(6, 13, 115, 7)),
              ICA.control = list(rseed = 1366))

# xi_Theta3
res9 <- bayes(fimfunc = ICAOD::FIM_sig_emax,
              lx = .001, ux = 500, k = 5, iter = 400,
              prior = uniform(c(4, 11, 100, 5), c(8, 15, 130, 9)),
              ICA.control = list(rseed = 1366))



# xi_Theta4
res10 <- bayes(fimfunc = ICAOD::FIM_sig_emax,
               lx = .001, ux = 500, k = 7, iter = 500,
               prior = uniform(c(4, 11, 100, 5), c(10, 18, 180, 11)),
               ICA.control = list(rseed = 13, ncount = 300, nimp = 30))

#######################################################################*
# Section 4.3 2-parameter Cox Proportional-Hazards Model for type one cenosred data ----
#######################################################################*
# The Fisher information matrix is available here with name FIM_2par_exp_censor1
# However, we should reparameterize the function to match the standard of the argument 'fimfunc'
myfim <- function(x, w, param)
  FIM_2par_exp_censor1(x = x, w = w, param = param, tcensor = 30)


res11 <- bayes(fimfunc = myfim, lx = 0, ux = 1, k = 2,
               iter = 80, prior = uniform(c(-3, -3), c(3, 3)),
               ICA.control = list(rseed = 1366))

res12 <- bayes(fimfunc = myfim, lx = 0, ux = 1, k = 4,
               iter = 150, prior = uniform(c(-11, -11), c(11, 11)),
               ICA.control = list(rseed = 1366))


#######################################################################*
# Section 4.3 2-parameter Cox Proportional-Hazards Model for type one cenosred data ----
#######################################################################*
# The Fisher information matrix is available here with name FIM_3par_exp_censor1
# However, we should reparameterize the function to match the standard of the argument 'fimfunc'
myfim2 <- function(x, w, param)
  FIM_3par_exp_censor1(x = x, w = w, param = param, tcensor = 30)


res13 <- bayes(fimfunc = myfim2, lx = 0, ux = 1, k = 3,
               iter = 100, prior = uniform(c(-3, -3, -3), c(3, 3, 3)),
               ICA.control = list(rseed = 1366))

res14 <- bayes(fimfunc = myfim2, lx = 0, ux = 1, k = 5,
               iter = 400, prior = uniform(c(-11, -11, -11), c(11, 11, 11)),
               ICA.control = list(rseed = 1366))



#################################*
# Section 4.4: DP-optimal designs ----
#################################*
# To show how to use formula interface, we use formula here.
p <- c(1, -2, 1, -1)
myprior <- uniform(p -1.5, p + 1.5)
myformula1 <- ~exp(b0+b1*x1+b2*x2+b3*x1*x2)/(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))
des_p0 <- bayescomp(formula = myformula1,
                    predvars = c("x1", "x2"),
                    parvars = c("b0", "b1", "b2", "b3"),
                    family = binomial(),
                    lx = c(-1, -1), ux = c(1, 1),
                    prior = myprior, iter = 50, k = 1,
                    prob = ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2)),
                    alpha = 0, ICA.control = list(rseed = 1366))



res15 <- bayescomp(formula = myformula1,
                   predvars = c("x1", "x2"),
                   parvars = c("b0", "b1", "b2", "b3"),
                   family = binomial(),
                   lx = c(-1, -1), ux = c(1, 1),
                   prior = myprior, iter = 1000, k = 7,
                   prob = ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2)),
                   alpha = .25, ICA.control = list(rseed = 1366, ncount = 40))


res16 <- bayescomp(formula = myformula1,
                   predvars = c("x1", "x2"),
                   parvars = c("b0", "b1", "b2", "b3"),
                   family = binomial(),
                   lx = c(-1, -1), ux = c(1, 1),
                   prior = myprior, iter = 1000, k = 7,
                   prob = ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2)),
                   alpha = .5, ICA.control = list(rseed = 13))

res17 <- bayescomp(formula = myformula1,
                   predvars = c("x1", "x2"),
                   parvars = c("b0", "b1", "b2", "b3"),
                   family = binomial(),
                   lx = c(-1, -1), ux = c(1, 1),
                   prior = myprior, iter = 1000, k = 7,
                   prob = ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2)),
                   alpha = .75, ICA.control = list(rseed = 1366,  ncount = 60))


res18 <- bayescomp(formula = myformula1,
                   predvars = c("x1", "x2"),
                   parvars = c("b0", "b1", "b2", "b3"),
                   family = binomial(),
                   lx = c(-1, -1), ux = c(1, 1),
                   prior = myprior, iter = 1000, k = 7,
                   prob = ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2)),
                   alpha = 1, ICA.control = list(rseed = 1366,  ncount = 60))


