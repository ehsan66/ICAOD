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

#############################################################*
# Exact  designs for the Compartmental model ----
##############################################################*
#Compare the exact design with the one reported in Gotwalt (2009)

# Gotwalt, C. M., B. A. Jones, and D. M. Steinberg (2009). Fast computation of designs robust to parameter uncertainty for nonlinear settings. Technometrics 51 (1), 88?95.

uni2  <- uniform(lower = c(.01884, .298), upper = c(.09884, 8.298))


exact <- bayes(formula = ~ b3 * (exp(-b1*x)- exp(-b2*x)),
               parvars = c("b1", "b2", "b3 = 21.8"),
               predvars = "x",
               lx = 0, ux = 48,
               k = 18,
               ICA.control = list(rseed = 1366, equal_weight = TRUE),
               iter = 500,
               prior = uni2, npar = 3)

exact$evol[[500]]$min_cost

beff(formula = ~ b3 * (exp(-b1*x)- exp(-b2*x)),
     parvars = c("b1", "b2", "b3 = 21.8"),
     predvars = "x", prior = uni2,
     x = c(rep(.2030, 5), rep(.9717, 4), rep(3.3439, 3), rep(20.1917, 6)),
     w = rep(1/18, 18),
     xopt = exact$evol[[500]]$x, wopt = exact$evol[[500]]$w)
# 0.9938858

