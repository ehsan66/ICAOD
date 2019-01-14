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


### This script calculates all the relative efficiencies reported in the Tables.

tun1 <- list(cubature = list(tol = 1e-5, maxEval = 50000))
#############################
# Section 4.1- 2PL model
############################
formula4.1 <- ~ 1/(1 + exp(-b *(x - a)))
predvars4.1 <- "x"
parvars4.1 <- c("a", "b")


des4.1 <- vector("list", 6)
des4.1[[1]]$x <- c(-3, -1.20829, 0, 1.20814, 3)
des4.1[[1]]$w <- c(.24701, .18305, .13988, .18309, .24702)
des4.1[[1]]$prior <- uniform(lower =  c(-3, .1), upper = c(3, 2))

des4.1[[2]]$x <- c(-2.41692, -1.16676, .04386, 1.18506, 2.40631)
des4.1[[2]]$w <- c(.26304, .18231, .14205, .16846, .24414)
des4.1[[2]]$prior <- student(mean =  c(0, 1), S   = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                             df = 3, lower =  c(-3, .1), upper = c(3, 2))

des4.1[[3]]$x <- c(-2.25540, -.76318, .54628, 2.16045)
des4.1[[3]]$w <- c(.31762, .18225, .18159, .31853)
des4.1[[3]]$prior <- normal(mu =  c(0, 1),
                            sigma = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                            lower =  c(-3, .1), upper = c(3, 2))

des4.1[[4]]$x <- c(-2.23013, -.66995, .67182, 2.23055)
des4.1[[4]]$w <- c(.31420, .18595, .18581, .31404)
des4.1[[4]]$prior <- normal(mu =  c(0, 1),
                            sigma = matrix(c(1, 0, 0, .5), nrow = 2),
                            lower =  c(-3, .1), upper = c(3, 2))

des4.1[[5]]$x <- c(-1.51175, .12043, 1.05272, 2.59691)
des4.1[[5]]$w <- c(.37679, .14078, .12676, .35567)
des4.1[[5]]$prior <- skewnormal(xi = c(0, 1),
                                Omega = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                                alpha = c(1, 0), lower =  c(-3, .1), upper = c(3, 2))


des4.1[[6]]$x <- c(-2.50914, -1.16780, -.36904, 1.29227)
des4.1[[6]]$w <- c(.35767, .11032, .15621, .37580)
des4.1[[6]]$prior <- skewnormal(xi = c(0, 1),
                                Omega = matrix(c(1, -0.17, -0.17, .5), nrow = 2),
                                alpha = c(-1, 0), lower =  c(-3, .1), upper = c(3, 2))


eff4.1 <- matrix(NA, 6, 6)
colnames(eff4.1) <- c("uni", "t", "norm1", "norm2", "skew1", "skew2")
rownames(eff4.1) <- colnames(eff4.1)
for (i in 1:6)
  for(j in 1:6)
    eff4.1[i, j] <- beff(formula = formula4.1,
                         predvars = predvars4.1,
                         parvars = parvars4.1,
                         family = binomial(),
                         prior = des4.1[[i]]$prior,
                         xopt = des4.1[[i]]$x,
                         wopt = des4.1[[i]]$w,
                         x = des4.1[[j]]$x,
                         w = des4.1[[j]]$w,
                         crt.bayes.control = tun1)
## Reported in Table 1 as eff_D
round(eff4.1[, 1], 3)
# 1.000 0.985 0.979 0.981 0.970 0.952
#############################
# Section 4.2- Emax model
############################
formula4.2 <- ~theta1 + (theta2 - theta1)*(x^theta4)/(x^theta4 + theta3^theta4)
predvars4.2 <- c("x")
parvars4.2 <- c("theta1", "theta2", "theta3", "theta4")

des4.2 <- vector("list", 4)


des4.2[[1]]$lb <- c(4, 11, 100, 5)
des4.2[[1]]$ub <- c(5, 12, 105, 6)
des4.2[[1]]$x <- c(0.00328, 84.75095, 123.83262, 500)
des4.2[[1]]$w <- c(.25, .25, .25, .25)
des4.2[[1]]$prior <- uniform(des4.2[[1]]$lb, des4.2[[1]]$ub)


###
des4.2[[2]]$lb <- c(4, 11, 100, 5)
des4.2[[2]]$ub <- c(6, 13, 115, 7)
des4.2[[2]]$x <- c(0.0882, 90.19404, 127.84294, 500)
des4.2[[2]]$w <-  c(.25, .25, .25, .25)
des4.2[[2]]$prior <- uniform(des4.2[[2]]$lb, des4.2[[2]]$ub)

#####
des4.2[[3]]$lb <- c(4, 11, 100, 5)
des4.2[[3]]$ub <- c(8, 15, 130, 9)
des4.2[[3]]$x <- c(0.36228, 94.59314, 113.66528, 138.3004, 500)
des4.2[[3]]$w <- c(.24322, .19415, .11571, .20332, .2436)
des4.2[[3]]$prior <- uniform(des4.2[[3]]$lb, des4.2[[3]]$ub)


des4.2[[4]]$lb <- c(4, 11, 100, 5)
des4.2[[4]]$ub <- c(10, 18, 180, 11)
des4.2[[4]]$x <- c(0.01849, 96.43274, 117.31918, 133.06849, 152.84389, 188.08142, 500)
des4.2[[4]]$w <- c(.20476, .11195, .1025, .09992, .13025, .13873, .21189)
des4.2[[4]]$prior <- uniform(des4.2[[4]]$lb, des4.2[[4]]$ub)


eff4.2 <- matrix(NA, 4, 4)
colnames(eff4.2) <- rownames(eff4.2) <- c("theta1", "theta2", "theta3", "theta4")

for (i in 1:4)
  for(j in 1:4)
    eff4.2[i, j] <- beff(formula = formula4.2,
                         predvars = predvars4.2,
                         parvars = parvars4.2,
                         prior = des4.2[[i]]$prior,
                         xopt = des4.2[[i]]$x,
                         wopt = des4.2[[i]]$w,
                         x = des4.2[[j]]$x,
                         w = des4.2[[j]]$w,
                         crt.bayes.control = tun1)
# Reported in Table 2 as eff_D
round(eff4.2[4, ], 3)
# 0.364 0.451 0.680 1.000

#############################
# Section 4.3- Censor model
############################
des4.3 <- vector("list", 4)
## two parameter
myfim1 <- function(x, w, param)
  FIM_2par_exp_censor1(x = x, w = w, param = param, tcensor = 30)
## three parameter
myfim2 <- function(x, w, param)
  FIM_3par_exp_censor1(x = x, w = w, param = param, tcensor = 1)


des4.3[[1]]$x <- c(0, 1)
des4.3[[1]]$w <- c(.5, .5)
des4.3[[1]]$prior <- uniform(c(-3, -3), c(3, 3))
des4.3[[1]]$fim <- myfim1
des4.3[[1]]$lb <- c(-3, -3)
des4.3[[1]]$ub <- c(3, 3)

des4.3[[2]]$x <- c(0,.30672, .55327, 1)
des4.3[[2]]$w <- c(.40548, .17300, .01427, .40724)
des4.3[[2]]$prior <- uniform(c(-11, -11), c(11, 11))
des4.3[[2]]$fim <- myfim1
des4.3[[2]]$lb <- c(-11, -11)
des4.3[[2]]$ub <- c(11, 11)

des4.3[[3]]$x <- c(0, .49042, 1)
des4.3[[3]]$w <- c(.33333, .33333, .33333)
des4.3[[3]]$prior <- uniform(c(-3, -3, -3), c(3, 3, 3))
des4.3[[3]]$fim <- myfim2
des4.3[[3]]$lb <- c(-3, -3, -3)
des4.3[[3]]$ub <- c(3, 3, 3)



des4.3[[4]]$x <- c(0, .18821, .45294, .6451, 1)
des4.3[[4]]$w <- c(.30267, .10198, .19726, .12909, .26900)
des4.3[[4]]$prior <- uniform(c(-11, -11, -11), c(11, 11, 11))
des4.3[[4]]$fim <- myfim2
des4.3[[4]]$lb <- c(-11, -11, -11)
des4.3[[4]]$ub <- c(11, 11, 11)


# Do not Run. The simulation takes more than one day!
len <- 20
eff4.3_locally <- matrix(NA, nrow = len^3, ncol = 4)

for(i in 1:4){
  b <- seq(des4.3[[i]]$lb[1], des4.3[[i]]$ub[1], length.out = len)
  b_two <- expand.grid(b, b)
  b_three <- expand.grid(b, b, b)
  if (length(des4.3[[i]]$lb) == 2){
    inipars <- b_two
  }else{
    inipars <- b_three
  }

  inipars <- as.matrix(inipars)
  for(j in 1:dim(inipars)[1]){
    iniparams <- inipars[j, ]
    temp <- locally(fimfunc = des4.3[[i]]$fim,
                    iter = 400,
                    k = length(iniparams), lx = 0, ux = 1, inipars = iniparams,
                    ICA.control = ICA.control(checkfreq = 200, plot_cost = FALSE, trace = FALSE))
    ## Efficiency of the locally optimal design with respect to the Bayesian optimal design under that prior
    eff4.3_locally[j, i] <- beff(fimfunc = des4.3[[i]]$fim,
                                 prior = des4.3[[i]]$prior,
                                 xopt = des4.3[[i]]$x,
                                 wopt = des4.3[[i]]$w,
                                 x = temp$evol[[400]]$x,
                                 w = temp$evol[[400]]$w,
                                 crt.bayes.control = tun1)
    temp <- NA
  }
  cat("\nsimulation number ", i, "is finished!, dim(b)[1] = ", dim(inipars)[1])
}

## reported in Table 3 as \bar{eff_D}
round(apply(eff4.3_locally, 2, mean, na.rm = TRUE), 3)
# 0.995 0.721 0.966 0.552

#############################
# Section 4.4- Compund criterion
############################
p <- c(1, -2, 1, -1)
prior4.4 <- uniform(p -1.5, p + 1.5)
formula4.4 <- ~exp(b0+b1*x1+b2*x2+b3*x1*x2)/(1+exp(b0+b1*x1+b2*x2+b3*x1*x2))
prob4.4 <- ~1-1/(1+exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))
predvars4.4 <-  c("x1", "x2")
parvars4.4 <- c("b0", "b1", "b2", "b3")
lb <- c(-1, -1)
ub <- c(1, 1)

des4.4 <- vector("list", 5)
des4.4[[1]]$x <- c(-1, 1)
des4.4[[1]]$w <- c(1)
des4.4[[1]]$alpha <- 0

des4.4[[2]]$x <- c(1, -.62534, .11405, -1, 1, .28175, -1, -1, 1, -1, -1, 1, 1, .09359)
des4.4[[2]]$w <- c(.08503, .43128, .01169, .14546, .05945, .08996, .17713)
des4.4[[2]]$alpha <- .25

des4.4[[3]]$x <- c(-1, .30193, 1, 1, .07411, -1, -.31952, -.08251, 1, -1, 1, -1, -1, 1)
des4.4[[3]]$w <- c(.09162, .10288, .15615, .13123, .01993, .22366, .27454)
des4.4[[3]]$alpha <- .5

des4.4[[4]]$x <- c(1, -1, .28274, 1, -1, -.19674, .03288, 1, -1, 1, -1, -.16751, 1, -1)
des4.4[[4]]$w <- c(.19040, .24015, .10011, .20527, .0388, .20075, .02452)
des4.4[[4]]$alpha <- .75

des4.4[[5]]$x <- c(1, -1, .26606, -.13370, 1, -.00887, -1, 1, -.2052, 1, 1, -1, -1, -1)
des4.4[[5]]$w <- c(.23020, .01612, .09546, .16197, .23675, .02701, .2325)
des4.4[[5]]$alpha <- 1

# reported in Table 4

# D-efficiency of the DP-optimal designs
# Reported in Table 4 as eff_D
# des4.4[[5]]$x and  des4.4[[5]]$w is the D-optimal design
beff(formula = formula4.4,
     predvars = predvars4.4,
     parvars = parvars4.4,
     family = binomial(),
     prior = prior4.4,
     xopt = des4.4[[5]]$x,
     wopt = des4.4[[5]]$w,
     x = des4.4[[2]]$x,
     w = des4.4[[2]]$w)
# 0.6443613

beff(formula = formula4.4,
     predvars = predvars4.4,
     parvars = parvars4.4,
     family = binomial(),
     prior = prior4.4,
     xopt = des4.4[[5]]$x,
     wopt = des4.4[[5]]$w,
     x = des4.4[[3]]$x,
     w = des4.4[[3]]$w)
# 0.9107611
beff(formula = formula4.4,
     predvars = predvars4.4,
     parvars = parvars4.4,
     family = binomial(),
     prior = prior4.4,
     xopt = des4.4[[5]]$x,
     wopt = des4.4[[5]]$w,
     x = des4.4[[4]]$x,
     w = des4.4[[4]]$w)
# 0.9875745

# must be one!
beff(formula = formula4.4,
     predvars = predvars4.4,
     parvars = parvars4.4,
     family = binomial(),
     prior = prior4.4,
     prob = prob4.4,
     type = "PA",
     xopt = des4.4[[5]]$x,
     wopt = des4.4[[5]]$w,
     x = des4.4[[5]]$x,
     w = des4.4[[5]]$w)
# 1
## P-efficiency
# reported in Table 4 as eff_P
# des4.4[[1]]$x and  des4.4[[1]]$w is the P-optimal design
beff(formula = formula4.4,
     predvars = predvars4.4,
     parvars = parvars4.4,
     family = binomial(),
     prior = prior4.4,
     prob = prob4.4,
     type = "PA",
     xopt = des4.4[[1]]$x,
     wopt = des4.4[[1]]$w,
     x = des4.4[[2]]$x,
     w = des4.4[[2]]$w)
# 0.8059524

beff(formula = formula4.4,
     predvars = predvars4.4,
     parvars = parvars4.4,
     family = binomial(),
     prior = prior4.4,
     prob = prob4.4,
     type = "PA",
     xopt = des4.4[[1]]$x,
     wopt = des4.4[[1]]$w,
     x = des4.4[[3]]$x,
     w = des4.4[[3]]$w)
# 0.6666876

beff(formula = formula4.4,
     predvars = predvars4.4,
     parvars = parvars4.4,
     family = binomial(),
     prior = prior4.4,
     prob = prob4.4,
     type = "PA",
     xopt = des4.4[[1]]$x,
     wopt = des4.4[[1]]$w,
     x = des4.4[[4]]$x,
     w = des4.4[[4]]$w)
# 0.5884567

beff(formula = formula4.4,
     predvars = predvars4.4,
     parvars = parvars4.4,
     family = binomial(),
     prior = prior4.4,
     prob = prob4.4,
     type = "PA",
     xopt = des4.4[[1]]$x,
     wopt = des4.4[[1]]$w,
     x = des4.4[[5]]$x,
     w = des4.4[[5]]$w)
#  0.5453367


