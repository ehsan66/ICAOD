# Locally D-optimal design for the 2PL model

# f <- function(x, a, b)
#   1/(1 + exp(-b * (x-a)))
# x <- seq(-3, 3, length.out = 100)
# y <- f(x, a = 0, b = .7)
# plot(x, y, type = "l")
setwd("../../my-papers/3-ICAOD JSS/")
library(ICAOD)
loc1 <- locally(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial",
                lx = -3, ux = 3, inipars = c(0, 1.25),
                iter = 100, k = 2,
                ICA.control = list(rseed = 1))


print(loc1)
loc1$evol[length(loc1$evol)]
plot(loc1)


pdf(file = "loc1.pdf")
plot(loc1)
dev.off()



#or
# loc1 <- locally(fimfunc = FIM_logistic,
#                   lx = -3, ux = 3, inipars = c(0, 1),  iter = 100, k = 2)

res1.2 <- locally(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                  parvars = c("a", "b"), family = binomial(),
                  lx = -3, ux = 3, inipars = c(0, 1.25), iter = 200, k = 2,
                  ICA.control = list(checkfreq = 50,
                                     stop_rule = "equivalence",
                                     stoptol = .99, rseed = 1))
print(res1.2)



loc3 <- locally(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = binomial(),
                lx = -3, ux = 3, inipars = c(0, 1.25), iter = 50,
                k = 3, x = c(0, 1.5, 2),
                ICA.control = list(rseed = 1, checkfreq = Inf))
print(loc3)


pdf(file = "loc3.pdf")
plot(loc3)
dev.off()

### minimax design for the 2PL model
min1 <- minimax(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial", iter = 100,
                lx = -3, ux = 3, lp = c(-3, .5), up = c(3, 2), k = 2,
                ICA.control = list(ncount = 40, rseed = 1,checkfreq = Inf,
                                   sym = TRUE, sym_point = 0),
                crt.minimax.control = list(optslist = list(maxeval = 200)))

print(min1)

pdf(file = "min1.pdf")
plot(min1)
dev.off()

## k = 3
ELB <- 0
k <- 3
time1 <- proc.time()
while(ELB < .99){
  min2 <- minimax(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                  parvars = c("a", "b"), family = "binomial",
                  #fimfunc = FIM_logistic,
                  lx = -3, ux = 3, lp = c(-3, .5), up = c(3, 2),
                  iter = 500, k = k,
                  ICA.control = list(ncount =40, rseed = 1, plot_cost = FALSE,
                                     sym = TRUE, sym_point = 0,
                                     checkfreq = Inf, trace = FALSE),
                  crt.minimax.control = list(optslist = list(maxeval = 200)))
  ELB <- min2$evol[[length(min2$evol)]]$sens$ELB
  # pdf(file = paste0("min2_",k, ".pdf"))
  # plot(min2)
  # dev.off()
  k <- k + 1
}
time2 <- proc.time() - time1

print(min2)

#1361/60

print(min2)
plot(min2)
pdf(file = "min2.1.pdf")
plot(min2[[1]])
dev.off()

# We use fimfunc argument instead of formula
# res1.4.1 <- minimax(fimfunc = FIM_logistic,
#                     lx = -3, ux = 3,
#                     lp = c(-3, .5), up = c(3, 2),
#                     iter = 500, k = 5,
#                     ICA.control = list(ncount =40, rseed = 1,
#                                        sym = TRUE, sym_point = 0),
#                     crt.minimax.control = list(optslist = list(maxeval = 200)))
#
# print(res1.4.1)
# ## CPU time comparison
# round((res1.4$arg$time[1] - res1.4.1$arg$time[1])/res1.4$arg$time[1], 2)


min3 <- minimax(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial",
                lx = -3, ux = 3, lp = c(-3, .5), up = c(3, 2),
                iter = 500, k = 5, n.grid = 11,
                ICA.control = list(rseed = 1, sym = TRUE,
                                   sym_point = 0, checkfreq = Inf))


### standardized maximin optimal designs for the 2pl model
# Also bring it to the minimax code
Dopt_2pl <- function(a, b){
  x <- c(a + (1/b) * 1.5434046, a - (1/b) * 1.5434046)
  return(list(x = x, w = c(.5, .5)))
}
min4 <- minimax(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial",
                standardized = TRUE, localdes = Dopt_2pl,
                lx = -3, ux = 3, lp = c(-3, .5), up = c(3, 2),
                iter = 500, k = 5,
                ICA.control = list(ncount =40, rseed = 1,
                                   sym = TRUE, sym_point = 0,
                                   checkfreq = Inf),
                crt.minimax.control = list(optslist = list(maxeval = 200)))



plot(min4)

# relative efficiency of the minimax design with respect to the standardized maximin design
# under standardized criterion
meff(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
     parvars = c("a", "b"), family = "binomial",
     lp = c(-3, .5), up = c(3, 2),
     x2 = min4$evol[[500]]$x, w2 = min4$evol[[500]]$w,
     x1 = min3$evol[[500]]$x, w1 = min3$evol[[500]]$w,
     standardized = TRUE,
     localdes = Dopt_2pl)


#### bayesian
bay1 <- bayes(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
              parvars = c("a", "b"), family = "binomial",
              lx = -3, ux = 3, iter = 600, k = 5,
              prior = uniform(lower = c(-3, .5), upper  = c(3, 2)),
              ICA.control = list(rseed = 1, checkfreq = Inf))
print(bay1)




pdf(file = "bay1.pdf")
plot(bay1)
dev.off()
# it is an example that shows how level affects the shape of the sensitivity function
#plot(res1.5, sens.bayes.control = list(method = "quadrature",
#                                       quadrature = list(level = 10)))


norm_prior <- normal(lower = c(-3, .5), upper  = c(3, 2), mu = c(0,1.25),
                     sigma = matrix(c(1, -.17, -.17, .5), 2))


bay2 <- bayes(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
              parvars = c("a", "b"), family = "binomial",
              lx = -3, ux = 3, iter = 1000, k = 4, prior = norm_prior,
              ICA.control = list(rseed = 1),
              crt.bayes.control = list(method = "quadrature",
                                       quadrature = list(level = 7)))


print(bay2)
plot(bay2, sens.bayes.control = list(method = "quadrature",
                                     quadrature = list(level = 7)))
pdf(file = "bay2.pdf")
plot(bay2, sens.bayes.control = list(method = "quadrature",
                                     quadrature = list(level = 7)))
dev.off()


beff(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
     parvars = c("a", "b"), family = "binomial",
     prior = norm_prior,
     x1 = bay1$evol[[600]]$x, w1 = bay1$evol[[600]]$w,
     x2 = bay2$evol[[600]]$x, w2 = bay2$evol[[600]]$w)



#### Optimum-on-Average or Robust designs
rob1 <- robust(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
               parvars = c("a", "b"), family = "binomial",
               lx = -3, ux = 3, iter = 1000, k = 4,
               prob = c(.25, .5, .25),
               parset = matrix(c(-2, 0, 2, 1.25, 1.25, 1.25), 3, 2),
               ICA.control = list(checkfreq = 50, stoptol = .999,
                                  stop_rule = "equivalence",
                                  rseed = 1))

print(rob1)

pdf(file = "rob1.pdf")
plot(rob1)
dev.off()

### exact
# uni_prior <- uniform(lower = c(0.01884, 0.298),
#                      upper = c(0.09884, 8.298))
# res1.8 <-  bayes(formula = ~ c*(exp( - a * t) - exp( - b * t)),
#                  predvars = "t", parvars = c("a", "b", "c = 21.8"),
#                  lx = 0, ux = 24, iter = 1000, k = 18, prior = uni_prior,
#                  ICA.control = list(rseed = 1, equal_weight = TRUE))
#
# #crt.bayes.control = list(method = "quadrature",
# #                          quadrature = list(level = 3)))
#
#
# print(res1.8)
#
#
#
# uni_prior2 <- uniform(lower = c(-3, 4, 5, -6, -2.5),
#                       upper = c(3, 10, 11, 0, 3.5))
# formula1 <- ~exp(b0 + b1 * x1 + b2 * x2 + b3 * x3 + b4 * x4)/(1+ exp(b0 + b1 * x1 + b2 * x2 + b3 * x3 + b4 * x4))
# res1.9 <-  bayes(formula = formula1,
#                  predvars = c("x1", "x2", "x3", "x4"),
#                  parvars = c("b0", "b1", "b2", "b3", "b4"),
#                  family = "binomial",
#                  lx = rep(-1, 4),  ux = rep(1, 4),
#                  iter = 400, k = 6,
#                  prior = uni_prior2,
#                  ICA.control = list(rseed = 1, equal_weight = TRUE),
#                  crt.bayes.control = list(method = "cubature",
#                                           cubature = list(maxEval = 5000)))
# #,
# #crt.bayes.control = list(method = "quadrature",
# #                           quadrature = list(level = 6)))
# res1.9$arg$time
#
#
#
# loc10 <-  bayes(formula = formula1,
#                   predvars = c("x1", "x2", "x3", "x4"),
#                   parvars = c("b0", "b1", "b2", "b3", "b4"),
#                   family = "binomial",
#                   lx = rep(-1, 4),  ux = rep(1, 4),
#                   iter = 1500, k = 10,
#                   prior = uni_prior2,
#                   ICA.control = list(rseed = 1, ncount = 120),
#                   crt.bayes.control = list(method = "quadrature",
#                                            quadrature = list(level = 4)))
# plot(loc10,
#      crt.bayes.control = list(method = "quadrature",
#                               quadrature = list(level = 3)))
# loc10$arg$time
# 10281/60

########
install.packages("VNM")
library(VNM)
###########################
## multiple
###########################
library(VNM)
Res1.MOPT1 <- MOPT(LB = log(0.001), UB = log(500),
                   P = c(60, 340, 107.14, 1),
                   lambda = c(1/3, 1/3), delta = 200, r = 30,
                   epsilon_w = 10^-7, verbose = TRUE)
summary(Res1.MOPT1)



# convert p1 to the parameters of the Hill model
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

mult1$evol[[100]]$x
summary(Res1.MOPT1)[1, ] # are the same

print(mult1)
pdf(file = "mult1.pdf")
plot(mult1)
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
# ###########################
# ## example 3
# ###########################
# Res2.MOPT <- MOPT(LB = log(0.001), UB = log(100),
#                   P = c(22, 16.8, 70, 1), lambda = c(1/3, 1/3), delta = 5, r = 30,
#                   verbose = TRUE)
#
# summary(Res2.MOPT)
#
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
# ###########################
# ## example 4
# ###########################
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
# ###########################
# ## example 5
# ###########################
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

###########################
## DP-optimal designs
###########################

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

#####
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

res4.1 <- locally(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                  parvars = c("a", "b"), family = "binomial",
                  lx = -1, ux = 1, inipars = c(0, 7),
                  iter = 100, k = 2,
                  crtfunc = c_opt, sensfunc = c_sens,
                  ICA.control = list(rseed = 1, checkfreq = Inf))
res4.1
pdf(file = "res4_1.pdf")
plot(res4.1)
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



res4.2 <- bayes(formula = ~1/(1 + exp(-b * (x-a))), predvars = "x",
                parvars = c("a", "b"), family = "binomial",
                lx = -1, ux = 1,
                prior = uniform(lower = c(-.3, 6), upper  = c(.3, 8)),
                iter = 100, k = 3,
                crtfunc = c_opt_vec,
                sensfunc = c_sens_vec,
                ICA.control = list(rseed = 1, ncount = 80),
                crt.bayes.control = list(cubature = list(maxEval = 50000)))

plot(res4.2,sens.bayes.control = list(method= "cubature",
                                      cubature = list(maxEval = 100)))

pdf(file = "res4_2.pdf")
plot(res4.2,  sens.bayes.control = list(method= "cubature",
                                        cubature = list(maxEval = 10000)))
dev.off()
optlist$algorithm


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
