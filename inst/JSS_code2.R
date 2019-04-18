

setwd("../../my-papers/3-ICAOD JSS/")
library(ICAOD)
############################################################################*
# 5.1: The Logisitc regression model with one single continuous variable ----
############################################################################*
log1 <- locally(formula = ~exp(b0 + b1 * x)/(1 + exp(b0 + b1 * x)),
                predvars = "x", parvars = c("b0", "b1"),
                family = "binomial", lx = 0, ux = 6, iter = 40, k = 2,
                inipars = c(-4, 1.3333), ICA.control = list(rseed = 1))


print(log1)
log1$evol[length(log1$evol)]
plot(log1)

pdf(file = "log1.pdf")
plot(log1)
dev.off()

leff(formula = ~exp(b0 + b1 * x)/(1 + exp(b0 + b1 * x)),
     predvars = "x", parvars = c("b0", "b1"),
     family = "binomial", inipars = c(-4, 1.3333),
     x1 = c(0:6), w1 = rep(1/7, 7),
     x2 = log1$evol[[40]]$x, w2 = log1$evol[[40]]$w)

round((1/0.7778719-1) * 100)

log2 <- locally(formula = ~exp(b0 + b1 * x)/(1 + exp(b0 + b1 * x)),
                predvars = "x", parvars = c("b0", "b1"),
                family = "binomial", lx = 0, ux = 6, iter = 40, k = 2,
                inipars = c(-4, 1.3333),
                ICA.control = list(rseed = 1,
                                   checkfreq = 20,
                                   stop_rule = "equivalence",
                                   stoptol = .99))



print(log2)




log3 <- locally(formula = ~exp(b0 + b1 * x)/(1 + exp(b0 + b1 * x)),
                predvars = "x", parvars = c("b0", "b1"),
                family = "binomial", lx = 0, ux = 6, iter = 40,
                x = c(1, 2, 3),
                inipars = c(-4, 1.3333),
                ICA.control = list(rseed = 1, checkfreq = Inf))
print(log3)
pdf(file = "log3.pdf")
plot(log3)
dev.off()


log4 <- minimax(formula = ~exp(b0 + b1 * x)/(1 + exp(b0 + b1 * x)),
                predvars = "x", parvars = c("b0", "b1"),
                family = "binomial",
                lx = 0, ux = 6, lp = c(-6, .5), up = c(-2, 2),
                iter = 200, k = 2,
                ICA.control = list(rseed = 1,
                                   checkfreq = 50,
                                   stop_rule = "equivalence",
                                   stoptol = .99),
                crt.minimax.control = list(optslist = list(maxeval = 200)))

print(log4)

pdf(file = "log4.pdf")
plot(log4)
dev.off()




log5 <- minimax(formula = ~exp(b0 + b1 * x)/(1 + exp(b0 + b1 * x)),
                predvars = "x", parvars = c("b0", "b1"),
                family = "binomial",
                lx = 0, ux = 6, lp = c(-6, .5), up = c(-2, 2),
                iter = 500, k = 3,
                ICA.control = list(rseed = 1,
                                   checkfreq = 50,
                                   stop_rule = "equivalence",
                                   stoptol = .99),
                crt.minimax.control = list(optslist = list(maxeval = 200)))

print(log5)

pdf(file = "log5.pdf")
plot(log5)
dev.off()


meff(formula = ~exp(b0 + b1 * x)/(1 + exp(b0 + b1 * x)),
     predvars = "x", parvars = c("b0", "b1"),
     family = "binomial",
     lp = c(-6, .5), up = c(-2, 2),
     x1 = c(0:6), w1 = rep(1/7, 7),
     x2 = log5$evol[[200]]$x, w2 = log5$evol[[200]]$w)
############################################################################*
# 5.2: The sigmiod Emax Model----
############################################################################*
### uniform prior
prior1 <- uniform(lower = c(4, 11, 100, 5), upper = c(8, 15, 130, 9))


sig1 <- bayes(formula = ~b1 + (b2-b1) * x ^b4/(x^b4 + b3^b4),
              predvars = "x",
              parvars = c("b1", "b2", "b3", "b4"),
              lx = .001, ux = 1000, k = 5, iter = 400,
              prior = prior1,
              ICA.control = list(rseed = 1, checkfreq = Inf))

print(sig1)
pdf(file = "sig1.pdf")
plot(sig1)
dev.off()



beff(formula = ~b1 + (b2-b1) * x ^b4/(x^b4 + b3^b4),
     predvars = "x",
     parvars = c("b1", "b2", "b3", "b4"),
     prior = prior1,
     x1 = c(.001,seq(100, 1000, by = 100)),
     w1 = rep(1/11, 11),
     x2 = sig1$evol[[400]]$x, w2 = sig1$evol[[400]]$w)
(1/0.3063289 - 1)*100

#### robust
parset1 <- matrix(c(4, 11, 100, 5,
                    5, 12, 110, 6,
                    6, 13, 120, 7,
                    8, 15, 130, 9,
                    12, 30, 160, 13),
                  nrow = 5, byrow = TRUE)

sig2 <- robust(formula = ~b1 + (b2-b1) * x ^b4/(x^b4 + b3^b4),
               predvars = "x",
               parvars = c("b1", "b2", "b3", "b4"),
               lx = .001, ux = 1000, k = 6, iter = 400,
               parset = parset1,
               prob = rep(1/5, 5),
               ICA.control = list(rseed = 1, checkfreq = Inf))
sig2
pdf(file = "sig2.pdf")
plot(sig2)
dev.off()
######################################################*
## an example to compare the hcubature algorithm
## with the traditional quadrature methods!
lp <- c(4, 11, 100, 5)
up <- c(8, 15, 500, 9)

prior2 <- normal(mu = (lp + up)/2,
                 sigma = diag(c(2, 2, 20^2, 3), 4),
                 lower = lp, upper = up)
# we adjust the tuning parameters because hcubature algorithm becomes
#  very slow for normal distributions.
#   it takes around 20 minutes on my laptop
sig3 <- bayes(formula = ~b1 + (b2-b1) * x^b4/(x^b4 + b3^b4),
              predvars = "x",
              parvars = c("b1", "b2", "b3", "b4"),
              lx = .001, ux = 1000, k = 5, iter = 400,
              prior = prior2,
              ICA.control = list(rseed = 1, checkfreq = Inf),
              crt.bayes.control = list(cubature = list(maxEval = 2000)),
              sens.bayes.control = list(cubature = list(maxEval = 2000)))

## We try finding the Bayesain Design with quadrature rule
sig4 <- bayes(formula = ~b1 + (b2-b1) * x^b4/(x^b4 + b3^b4),
              predvars = "x",
              parvars = c("b1", "b2", "b3", "b4"),
              lx = .001, ux = 1000, k = 5, iter = 400,
              prior = prior2,
              ICA.control = list(rseed = 1, checkfreq = Inf),
              crt.bayes.control = list(method = "quadrature",
                                       quadrature = list(level = 10,
                                                         type = "GHe")),
              sens.bayes.control = list(method = "cubature",
                                        cubature = list(maxEval = 5000)))
# Even a quadrature rule with 10 points will not produce
#   the optimal design after 90 min on my system!!
# This is an example to show why the authors prefer the hcubature algorithm
#  over the traditional quadrature methods!
######################################################*



############################################################################*
# 5.3: c-optimal designs----
############################################################################*
c_opt <-function(x, w, a, b, fimfunc){
  gam <- log(.95/(1-.95))
  M <- fimfunc(x = x, w = w, a = a, b = b)
  c <- matrix(c(1, -gam * b^(-2)), nrow = 1)
  B <- t(c) %*% c
  sum(diag(B %*% solve(M)))
}

c_sens <- function(xi_x, x, w, a, b, fimfunc){
  gam <- log(.95/(1-.95))
  M <- fimfunc(x = x, w = w, a = a, b = b)
  M_inv <- solve(M)
  M_x <- fimfunc(x = xi_x, w = 1, a = a, b = b)
  c <- matrix(c(1, -gam * b^(-2)), nrow = 1)
  B <- t(c) %*% c
  sum(diag(B %*% M_inv %*% M_x %*%  M_inv)) - sum(diag(B %*% M_inv))
}

twoPL1 <- locally(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                  parvars = c("a", "b"), family = "binomial",
                  lx = -1, ux = 1, inipars = c(0, 7),
                  iter = 100, k = 2,
                  crtfunc = c_opt, sensfunc = c_sens,
                  ICA.control = list(rseed = 1, checkfreq = Inf))
twoPL1
pdf(file = "twoPL1.pdf")
plot(twoPL1)
dev.off()

## bayesian
c_opt_vec <-function(x, w, a, b, fimfunc){
  gam <- log(.95/(1-.95))
  M <- fimfunc(x = x, w = w, a = a, b = b)
  B <- sapply(1:length(M), FUN = function(i)
    matrix(c(1, -gam * b[i]^(-2)), ncol= 1) %*%
      matrix(c(1, -gam * b[i]^(-2)), nrow = 1), simplify = FALSE)
  sapply(1:length(M), FUN = function(i)
    sum(diag(B[[i]] %*% solve(M[[i]]))))
}

c_sens_vec <- function(xi_x, x, w, a, b, fimfunc){
  gam <- log(.95/(1-.95)) # LD .95
  M <- fimfunc(x = x, w = w, a = a, b = b)
  M_inv <- lapply(M , FUN = function(FIM) solve(FIM))
  M_x <- fimfunc(x = xi_x, w = 1, a = a, b = b)
  B <- sapply(1:length(M), FUN = function(i)
    matrix(c(1, -gam * b[i]^(-2)), ncol= 1) %*%
      matrix(c(1, -gam * b[i]^(-2)), nrow = 1), simplify = FALSE)
  sapply(1:length(M), FUN = function(i)
    sum(diag(B[[i]] %*% M_inv[[i]] %*% M_x[[i]] %*% M_inv[[i]])) -
      sum(diag(B[[i]] %*% M_inv[[i]])))
}



twoPL2 <- bayes(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial",
                lx = -1, ux = 1,
                prior = uniform(lower = c(-.3, 6), upper  = c(.3, 8)),
                iter = 100, k = 3,
                crtfunc = c_opt_vec,
                sensfunc = c_sens_vec,
                ICA.control = list(rseed = 1, ncount = 60, checkfreq = Inf),
                sens.bayes.control = list(cubature = list(maxEval = 100)))


print(twoPL2)

pdf(file = "twoPL2.pdf")
plot(twoPL2,  sens.bayes.control = list(method= "cubature",
                                        cubature = list(maxEval = 10000)))
dev.off()



# A_bayes <- bayes(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
#                  parvars = c("a", "b"), family = "binomial",
#                  lx = -1.2, ux = 1.2,
#                  prior = uniform(lower = c(-.3, 6.9), upper  = c(.3, 7.1)),
#                  iter = 200, k = 3,
#                  crtfunc = c_opt_vec,
#                  sensfunc = c_sens_vec,
#                  ICA.control = list(rseed = 1, ncount = 100),
#                  crt.bayes.control = list(method= "cubature",
#                                           cubature = list(maxEval = 10000)))


###################################################*
# Supplementary material. Section 1
# Multiple-Objective Optimal Designs for the 4-Parameter Hill model----
###################################################*
install.packages("VNM")
library(VNM)

# convert p1 to the parameters of the Hill model
# the following parameters are used in the example of the VNM package
# we should convert them to a, b, c, d, the parameters of the Hill model
p1 <- c(60, 340, 107.14, 1)
theta1 <- c(p1[2], -p1[4],NA , p1[1])
theta1[3] = -log(p1[3]) * theta1[2]
c <- theta1[4]
b <- -theta1[2]
d <- theta1[1] + c
a <- exp(theta1[3]/b)

mult1 <- multiple(minDose = .001, maxDose = 500,
                  inipars = c(107.14, 1, 60, 400), k = 4, iter = 500,
                  lambda = c(1/3, 1/3, 1/3), delta = 200,
                  ICA.control = list(rseed = 1, ncount = 40,
                                     stop_rule = "equivalence",
                                     checkfreq = 100, stoptol = .99))


print(mult1)
pdf(file = "mult1.pdf")
plot(mult1)
dev.off()



library(VNM)
t1 <- proc.time()
Res1.MOPT1 <- MOPT(LB = log(0.001), UB = log(500),
                   P = c(60, 340, 107.14, 1),
                   lambda = c(1/3, 1/3), delta = 200, r = 30,
                   epsilon_w = 10^-7, verbose = TRUE)
t2 <- proc.time()-t1
summary(Res1.MOPT1)
plot(Res1.MOPT1)
mult1$evol[[100]]$x
summary(Res1.MOPT1)[1, ] # are the same

pdf(file = "Res1_MOPT1.pdf")
plot(Res1.MOPT1)
dev.off()


mult2 <- multiple(minDose = .001, maxDose = 500, iter = 100,
                  inipars = c(107.14, 1, 60, 400), k = 4,
                  lambda = c(1/3, 1/3, 1/3), delta = 200,
                  x = c(.001, 20, 200, 500),
                  ICA.control = list(rseed = 1, checkfreq = Inf))
print(mult2)
plot(mult2)


pdf(file = "mult2.pdf")
plot(mult2)
dev.off()

Res2.MOPT <- MOPT(LB = log(0.001), UB = log(100),
                  P = c(22, 16.8, 70, 1), lambda = c(1/3, 1/3), delta = 5, r = 30,
                  verbose = TRUE)

summary(Res2.MOPT)

#
# p2 <-  c(22, 16.8, 70, 1)
# theta2 <- c(p2[2], -p2[4],NA , p2[1])
# theta2[3] = -log(p2[3]) * theta2[2]
#
# res9 <- multiple(minDose = log(.001), maxDose = log(100),
#                  inipars = theta2, k = 4, lambda = c(1/3, 1/3, 1/3),
#                  delta = 5,
#                  iter = 500,
#                  ICA.control = list(rseed = 1366, ncount = 100,
#                                     stop_rule = "equivalence",
#                                     checkfreq = 100, stoptol = .99))
#
# ###########################*
# ## example 4
# ###########################*
# res10 <-  multiple(minDose = log(.001),
#                    maxDose = log(500),
#                    inipars = theta1, k = 4,
#                    lambda = c(1/3, 1/3, 1/3),
#                    delta = 200,
#                    iter = 500,
#                    x = res8$evol[[length(res8$evol)]]$x,
#                    ICA.control = list(rseed = 1366, ncount = 100,
#                                       stop_rule = "equivalence",
#                                       checkfreq = 100, stoptol = .99))
#
# ###########################*
# ## example 5
# ###########################*
# res11 <-  multiple(minDose = log(.001),
#                    maxDose = log(100),
#                    inipars = theta2, k = 4,
#                    lambda = c(1/3, 1/3, 1/3),
#                    delta = 5,
#                    iter = 500,
#                    x = res9$evol[[length(res9$evol)]]$x,
#                    ICA.control = list(rseed = 1366, ncount = 100,
#                                       stop_rule = "equivalence",
#                                       checkfreq = 100, stoptol = .99))
#
#

###################################################*
# Supplementary material. Section 2
# DP-optimal designs----
###################################################*

p0 <- c(1, -2, 1, -1)
form3 <- ~exp(b0+b1*x1+b2*x2+b3*x1*x2)/(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))
prob3 <- ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))
pred3 <-  c("x1", "x2")
par3 <- c("b0", "b1", "b2", "b3")



# set checkfreq = Inf to ask for equivalence theorem at final step.
dp <- locallycomp(formula = form3, predvars = pred3, parvars = par3,
                  family = binomial(), prob = prob3, alpha = .5,
                  lx = c(-1, -1), ux = c(1, 1), k = 4, inipars = p0,
                  iter = 300, ICA.control = ICA.control(rseed = 1, checkfreq = Inf))
print(dp)

# or plot the sensitivity function with rgl
plot(res3.1, plot_3d = "rgl")

# or plot the sensitivity function with rgl
pdf(file = "dp_%d.pdf", onefile=FALSE)
plot(dp)
dev.off()
#
# res3.2 <- bayescomp(formula = form3, predvars = pred3, parvars = par3,
#                     family = binomial(), prob = prob3,
#                     lx = c(-1, -1), ux = c(1, 1), alpha = .5,
#                     k = 4, prior = uniform(p0 -1.5, p0 + 1.5),
#                     iter = 300, ICA.control = ICA.control(checkfreq = Inf),
#                     sens.bayes.control = list(cubature = list(maxEval = 1000)),
#                     crt.bayes.control = list(cubature = list(maxEval = 1000)))
# print(res3.1)
#


