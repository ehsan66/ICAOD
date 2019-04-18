## This script verifies Bayesian optimality of all the presented designs in the submitted paper.

## to install and load the package ICAOD
install_and_load <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

install_and_load("ICAOD")
install_and_load("rgl")

## the tuning parameters for the hcubature algorithm
# I decreased the maxEval to have faster code
tun1 <- list(cubature = list(tol = 1e-5, maxEval = 10000))

#############################**
# Section 4.1 Table 1 ----
############################**
formula4.1 <- ~ 1/(1 + exp(-b *(x - a)))
predvars4.1 <- "x"
parvars4.1 <- c("a", "b")

##### xi_uni
uni <- uniform(lower =  c(-3, .1), upper = c(3, 2))

res4.1_ex1 <- sensbayes(formula = formula4.1, predvars = predvars4.1, parvars = parvars4.1,
                        family = binomial(), lx = -3, ux = 3,
                        prior = uni,
                        x = c(-3, -1.20829, 0, 1.20814, 3),
                        w = c(.24701, .18305, .13988, .18309, .24702),
                        sens.bayes.control = tun1)
# 11 seconds

##### xi_t
stud <- student(mean =  c(0, 1), S   = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                df = 3, lower =  c(-3, .1), upper = c(3, 2))

res4.1_ex2 <- sensbayes(formula = formula4.1, predvars = predvars4.1, parvars = parvars4.1,
                        family = binomial(), lx = -3, ux = 3,
                        prior = stud,
                        x = c(-2.41692, -1.16676, .04386, 1.18506, 2.40631),
                        w = c(.26304, .18231, .14205, .16846, .24414),
                        sens.bayes.control = tun1)
# 27 seconds

##### xi_norm1
norm1 <- normal(mu =  c(0, 1),
                sigma = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                lower =  c(-3, .1), upper = c(3, 2))

res4.1_ex3 <- sensbayes(formula = formula4.1, predvars = predvars4.1, parvars = parvars4.1,
                        family = binomial(), lx = -3, ux = 3,
                        prior = norm1,
                        x = c(-2.25540, -.76318, .54628, 2.16045),
                        w = c(.31762, .18225, .18159, .31853),
                        sens.bayes.control = tun1)
# 21 seconds

##### xi_norm2
norm2 <- normal(mu =  c(0, 1),
                sigma = matrix(c(1, 0, 0, .5), nrow = 2),
                lower =  c(-3, .1), upper = c(3, 2))

res4.1_ex4 <- sensbayes(formula = formula4.1, predvars = predvars4.1, parvars = parvars4.1,
                        family = binomial(), lx = -3, ux = 3,
                        prior = norm2,
                        x = c(-2.23013, -.66995, .67182, 2.23055),
                        w = c(.31420, .18595, .18581, .31404),
                        sens.bayes.control = tun1)
# 16 seconds

##### xi_skew1
skew1 <- skewnormal(xi = c(0, 1),
                    Omega = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                    alpha = c(1, 0), lower =  c(-3, .1), upper = c(3, 2))

res4.1_ex5 <- sensbayes(formula = formula4.1, predvars = predvars4.1, parvars = parvars4.1,
                        family = binomial(), lx = -3, ux = 3,
                        prior = skew1,
                        x = c(-1.51175, .12043, 1.05272, 2.59691),
                        w = c(.37679, .14078, .12676, .35567),
                        sens.bayes.control = tun1)
# 23 seconds

##### xi_skew2
skew2 <- skewnormal(xi = c(0, 1),
                    Omega = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                    alpha = c(-1, 0), lower =  c(-3, .1), upper = c(3, 2))

res4.1_ex6 <- sensbayes(formula = formula4.1, predvars = predvars4.1, parvars = parvars4.1,
                        family = binomial(), lx = -3, ux = 3,
                        prior = skew2,
                        x = c(-2.50914, -1.16780, -.36904, 1.29227),
                        w = c(.35767, .11032, .15621, .37580),
                        sens.bayes.control = tun1)
# 16 seconds
#############################*
# Section 4.2 Table 2 ----
############################*
formula4.2 <- ~theta1 + (theta2 - theta1)*(x^theta4)/(x^theta4 + theta3^theta4)
predvars4.2 <- c("x")
parvars4.2 <- c("theta1", "theta2", "theta3", "theta4")

##### 4.2 Example 1
# Theta1
lb4.2_ex1 <- c(4, 11, 100, 5)
ub4.2_ex1 <- c(5, 12, 105, 6)

res4.2_ex1 <- sensbayes(formula = formula4.2,
                        predvars = predvars4.2, parvars = parvars4.2,
                        x = c(0.00328, 84.75095, 123.83262, 500),
                        w = c(.25, .25, .25, .25),
                        prior = uniform(lb4.2_ex1, ub4.2_ex1),
                        lx = .001, ux = 500,
                        sens.bayes.control = tun1)
# 9 seconds

##### 4.2 Example 2
# Theta2
lb4.2_ex2 <- c(4, 11, 100, 5)
ub4.2_ex2 <- c(6, 13, 115, 7)

res4.2_ex2 <- sensbayes(formula = formula4.2,
                        predvars = predvars4.2, parvars = parvars4.2,
                        x = c(0.0882, 90.19404, 127.84294, 500),
                        w = c(.25, .25, .25, .25),
                        prior = uniform(lb4.2_ex2, ub4.2_ex2),
                        lx = .001, ux = 500,
                        sens.bayes.control = tun1)
# 10 seconds

##### 4.2 Example 3
# Theta 3
lb4.2_ex3 <- c(4, 11, 100, 5)
ub4.2_ex3 <- c(8, 15, 130, 9)


res4.2_ex3 <- sensbayes(formula = formula4.2,
                        predvars = predvars4.2, parvars = parvars4.2,
                        x = c(0.36228, 94.59314, 113.66528, 138.3004, 500),
                        w = c(.24322, .19415, .11571, .20332, .2436),
                        prior = uniform(lb4.2_ex3, ub4.2_ex3),
                        lx = .001, ux = 500,
                        sens.bayes.control = tun1)
# 25 seconds

##### 4.2 Example 4
# Theta4
lb4.2_ex4 <- c(4, 11, 100, 5)
ub4.2_ex4 <- c(10, 18, 180, 11)

res4.2_ex4_cub <- sensbayes(formula = formula4.2,
                        predvars = predvars4.2, parvars = parvars4.2,
                        x = c(0.01849, 96.43274, 117.31918, 133.06849, 152.84389, 188.08142, 500),
                        w = c(.20476, .11195, .1025, .09992, .13025, .13873, .21189),
                        prior = uniform(lb4.2_ex4, ub4.2_ex4), lx = .001, ux = 500,
                        sens.bayes.control = tun1)
# around 1 min

#############################*
# Section 4.3 Table 3 ----
############################*
#### Cox proportional-hazard model

# two paramter
myfim1 <- function(x, w, param)
  FIM_2par_exp_censor1(x = x, w = w, param = param, tcensor = 30)

res4.3_2par_d3 <- sensbayes(fimfunc = myfim1, lx = 0, ux = 1,
                            prior = uniform(c(-3, -3), c(3, 3)),
                            x = c(0, 1), w = c(.5, .5),
                            sens.bayes.control = tun1)
# 7 seconds

res4.3_2par_d11 <- sensbayes(fimfunc = myfim1, lx = 0, ux = 1,
                             prior = uniform(c(-11, -11), c(11, 11)),
                             x = c(0,.30672, .55327, 1),
                             w = c(.40548, .17300, .01427, .40724),
                             sens.bayes.control = list(cubature = list(tol = 1e-6, maxEval = 100000)))
# 180s

# three paramter
myfim2 <- function(x, w, param)
  FIM_3par_exp_censor1(x = x, w = w, param = param, tcensor = 30)

res4.3_3par_d3 <- sensbayes(fimfunc = myfim2,
                            lx = 0, ux = 1,
                            prior = uniform(c(-3, -3, -3), c(3, 3, 3)),
                            x = c(0, .49042, 1), w = c(.333, .333, .333),
                            sens.bayes.control = tun1)
# 14 seconds

res4.3_3par_d11 <- sensbayes(fimfunc = myfim2, lx = 0, ux = 1,
                             prior = uniform(c(-11, -11, -11),
                                             c(11,11, 11)),
                             x = c(0, .18821, .45294, .6451, 1),
                             w = c(.30267, .10198, .19726, .12909, .26900),
                             sens.bayes.control = tun1)
# around 973 seconds!
# ELB 0.9986973
# Let us try it with less conservative tuning parameters
res4.3_3par_d11 <- sensbayes(fimfunc = myfim2, lx = 0, ux = 1,
                             prior = uniform(c(-11, -11, -11),
                                             c(11,11, 11)),
                             x = c(0, .18821, .45294, .6451, 1),
                             w = c(.30267, .10198, .19726, .12909, .26900),
                             sens.bayes.control = list(cubature = list(tol = 1e-4, maxEval = 10000)))
# 125.223 seconds!
# ELB 0.9982586
# So it really does not make much difference. You can use both sets of tuning parameters
# This is only true for the equivalence theorem
############################*
# Section 4.4 Table 4 ----
############################*
p <- c(1, -2, 1, -1)
prior4.4 <- uniform(p -1.5, p + 1.5)
formula4.4 <- ~exp(b0+b1*x1+b2*x2+b3*x1*x2)/(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))
prob4.4 <- ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))
predvars4.4 <-  c("x1", "x2")
parvars4.4 <- c("b0", "b1", "b2", "b3")

## we reduced the tolerance here to plot the derivative functions in less CPU time.
## But tol = 1e-5 is also applicable here, with higher CPU time.
tun2 <- list(cubature = list(tol = 1e-3, maxEval = 50000))
# In tthe paper, we used 1e-5 to be very sure of the optimality confirmation.
# in this case, it can take 15 min.
# with tol=1e-5, the ELB is approximately 1

### alpha= .25
# 1 min on my system when tol = 1e-3
# use package rgl
res4.4_alpha.25 <- sensbayescomp(formula = formula4.4,
                                 predvars = predvars4.4, parvars = parvars4.4,
                                 family = binomial(),
                                 lx = c(-1, -1), ux = c(1, 1),
                                 x = c(1, -.62534, .11405, -1, 1, .28175, -1,
                                       -1, 1, -1, -1, 1, 1, .09359),
                                 w  = c(.08503, .43128, .01169, .14546, .05945, .08996, .17713),
                                 prior = prior4.4,
                                 prob = prob4.4,
                                 alpha = .25,
                                 sens.bayes.control = tun2, plot_3d = "rgl")


### alpha= .5
# 30 seconds on my system when tol = 1e-3
res4.4_alpha.5 <- sensbayescomp(formula = formula4.4,
                                predvars = predvars4.4, parvars = parvars4.4,
                                family = binomial(),
                                lx = c(-1, -1), ux = c(1, 1),
                                x = c(-1, .30193, 1, 1, .07411, -1, -.31952,
                                      -.08251, 1, -1, 1, -1, -1, 1),
                                w  = c(.09162, .10288, .15615, .13123, .01993, .22366, .27454),
                                prior = prior4.4,
                                prob = prob4.4,
                                alpha = .5, plot_3d = "rgl",
                                sens.bayes.control = tun2)

### alpha = .75
# 40 seconds on my system when tol = 1e-3
res4.4_alpha.75 <- sensbayescomp(formula = formula4.4,
                                 predvars = predvars4.4, parvars = parvars4.4,
                                 family = binomial(),
                                 lx = c(-1, -1), ux = c(1, 1),
                                 x = c(1, -1, .28274, 1, -1, -.19674, .03288,
                                       1, -1, 1, -1, -.16751, 1, -1),
                                 w  = c(.19040, .24015, .10011, .20527, .0388, .20075, .02452),
                                 prior = prior4.4,
                                 prob = prob4.4,
                                 alpha = .75,
                                 sens.bayes.control = tun2)

### alpha = 1
# 30 seconds on my system when tol = 1e-3
res4.4_alpha.1 <- sensbayescomp(formula = formula4.4,
                                predvars = predvars4.4, parvars = parvars4.4,
                                family = binomial(),
                                lx = c(-1, -1), ux = c(1, 1),
                                x = c(1, -1, .26606, -.13370, 1, -.00887, -1,
                                      1, -.2052, 1, 1, -1, -1, -1),
                                w  = c(.23020, .01612, .09546, .16197, .23675, .02701, .2325),
                                prior = prior4.4,
                                prob = prob4.4,
                                alpha = 1,
                                sens.bayes.control = tun2)
