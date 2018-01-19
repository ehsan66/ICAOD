devtools::load_all(pkg = "D:\\Dropbox\\PhD Projects\\Minimax optimal design\\dod")

#######################################################################################################################
########## for equivalence

#############################################################
## check locally optimality: lp = up and type = "locally"
inipar <- c(2, 3)
equivalence (fimfunc = "FIM_logistic",
             x = c(1.485526, 2.51446 ),
             w = c(.5, .5),
             lx = -5, ux = 5,
             lp = inipar, up = inipar,
             type = "locally")

test <- mica(fimfunc = "exp_2par", lx = 0, ux = 1,
             lp = c(1, 1), up = c(1, 10), iter = 400,
             k = 3, type = "standardized",
             control = list(seed = 215))



## there is no analytical solution for locally optimal design for this model
## gosolnp automatically will be used to find the locally optimal design in the denominator
## of standardized criterion. Becasue, it is two-level nested optimization (first level on parameter space)
## ans second level on design space) it takes so long to find 'all_optima' and construct 'answerign' set.
equivalence (fimfunc = "FIM_power_logistic",
             x = c(-4.5515, 0.2130, 2.8075),
             w = c(0.4100, 0.3723, 0.2177),
             lx = -5, ux = 5,
             lp = c(0, 1), up = c(3, 1.5),
             type = "standardized",
             s = .2)


## standardized maximin D-optimal design does not depend on theta0 and theta1, so we fix them
## locally D-optimal design has a closed-form which is defined internally
equivalence (fimfunc = "FIM_loglin",
             x = c(0, 4.2494, 17.0324, 149.9090),
             w = c(0.3204, 0.1207, 0.2293, 0.3296),
             lx = 0, ux = 150,
             lp = c(2, 2, 1), up = c(2, 2, 15),
             type = "standardized")


### when a design point is of two dimension
# using package rgl (rgl must be installed for plot)
equivalence (fimfunc = "FIM_mixed_inhibition",
             x = c(3.4614, 4.2801, 30, 30, 0, 3.1426, 0, 4.0373 ),
             w = rep(1/4, 4),
             lx = c(0, 0), ux = c(30, 60),
             lp = c(7, 4, 2, 4), up = c(7, 5, 3, 5),
             type = "standardized", plot_3d = "rgl")
# here the design points are x1 = c(3.4614, 0), x2 = c(4.2801, 3.1426), x3 = c(30, 0), x4 = c(30, 4.0373)


equivalence (fimfunc = "FIM_comp_inhibition",
             x = c(3.4432, 30, 30, 0, 0,  18.8954),
             w = rep(1/3, 3),
             lx = c(0, 0), ux = c(30, 60),
             lp = c(7, 4,  4), up = c(7, 5,  5),
             type = "standardized")


locally_det<- function(param,  auxiliary){
  ## param is the vector of theta = (theta0, theta1, theta2)
  ux <- 0
  lx <- 150
  xstar <- (ux + param[3]) * (lx + param[3]) * (log(ux + param[3]) - log(lx + param[3]))/(ux - lx) - param[3]
  denominator <- det(FIM_loglin(x = c(lx, xstar, ux) , w = rep(1/3, 3), param = param))
  return(denominator)
}

equivalence (fimfunc = "FIM_loglin",
             x = c(0, 4.2494, 17.0324, 149.9090),
             w = c(0.3204, 0.1207, 0.2293, 0.3296),
             lx = 0, ux = 150,
             lp = c(2, 2, 1), up = c(2, 2, 15),
             locally = locally_det,
             type = "standardized")

#########################################################################################################################
#######################################################################################
#######################################################################################
######## optim on the average

test <- on_average_ica (fimfunc = "FIM_logistic", lx = -5, ux = 5, prior = rep(1/4, 4),
                        param = matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2),
                        iter = 200, k =3)

plot(test)
print(test)


test <- on_average_ica (fimfunc = "FIM_logistic", lx = -5, ux = 5, prior = rep(1/4, 4),
                        param = matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2),
                        iter = 200, k =3,
                        control = list(stop_rule = "equivalence",
                                       stoptol = .9995, equivalence_every = 100))

equivalence_ave(fimfunc ="FIM_logistic",lx = -5, ux = 5, x = c(0.2603688, 1, 1.739631),
                       w = c(0.2750147, 0.4499705, 0.2750148),  prior = c(.25, .25, .25, .25),
                       param =  matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2))







#####################################################################################################################
### for multiobjective optimal design

## An example how to create the design in Hyun and Wong (2015)
## An initial guess from Table 1:
Theta1 <- c(1.563, 1.790, 8.442, 0.137)

#########################################################
## Table 2 first row
# creating optimal design for estimating parameters: \xi_D
test <- multica_4pl(lx = log(.001),
                    ux = log(1000),
                    param = Theta1,
                    k = 4,
                    lambda = c(1, 0, 0),
                    delta = -1,
                    iter = 150,
                    control = list(seed = 1366, plot_cost = TRUE))
plot(test)
########################################################
# creating optimal design for estimating ED50: \xi_{ED50}
test2 <- multica_4pl(lx = log(.001),
                     ux = log(1000),
                     param = Theta1,
                     k = 3,
                     lambda = c(0, 1, 0),
                     delta = -1,
                     iter = 700,
                     control = list(seed = 1366, nimp = 20, ncount = 200, equivalence_every = 500))


########################################################
# creating optimal design for estimating MED: \xi_{MED}
test  <- multica_4pl(lx = log(.001),
                     ux = log(1000),
                     param = Theta1,
                     k = 3,
                     lambda = c(0, 0, 1),
                     delta = -1,
                     iter = 3000,
                     control = list(seed = 1366, nimp = 50, ncount = 500))



#######################################################
## finding multiple objective optimal design: example 1, Table 3
res1 <- multica_4pl(lx = log(.001),
                    ux = log(1000),
                    param = Theta1,
                    k = 4, lambda = c(0.05, 0.05, .90),
                    delta = -1, iter = 400,
                    control = list(seed = 1366))

plot(res1)

#######################################################
## finding multiple objective optimal design: example 2, Table 3
res2 <- multica_4pl(lx = log(.001),
                    ux = log(1000),
                    param = c(16.8, -1, 4.248, 22),
                    k = 4, lambda = c(0, 0.1, .9),
                    delta = 5, iter = 200,
                    control = list(seed = 1366))
plot(res2)
#######################################################################################
######################################################################################


#######################################################################################
## how to transfer from Hill model to 4-parameter logistic model
## parameters for Hill model
a <- 0.008949  # ED50
b <- -1.79 # Hill constant
c <- 0.137 # lower limit
d <- 1.7 # upper limit
D <- c(.001, 1000) ## dose in mg
## Hill_para is c(a, b, c, d)


res2 <- multica_4pl(lx = log(.001),
                    ux = log(1000),
                    param =  c(d - c, -b, b * log(a), c),
                    k = 4, lambda = c(0.05, 0.05, .90),
                    delta = -1, iter = 400,
                    control = list(seed = 1366))
exp(res2$evol[[length(res2$evol)]]$x) # dose level in mg
#######################################################################################



##########################################################################################
## equivalence_multiple examples

## verfying the design in Table 2 of Hyun and Wong (2015), first row, fisrt column.
Theta1 <- c(1.563, 1.790, 8.442, 0.137)
equivalence_multiple (x = c(log(.001), -5.21, -4.08, log(1000)),
                      w = c(.25, .25, .25, .25),
                      lx = log(.001), ux = log(1000),
                      param = Theta1,
                      lambda = c(1, 0, 0),
                      delta = -1)


#########################################################################################
## examples fof using this function for c-optimal designs
# first row second column: c-optimal design for estimating ED50
equivalence_multiple (x = c(log(.001), -4.80, log(1000)),
                      w = c(.276, .500, .224),
                      lx = log(.001), ux = log(1000),
                      param = Theta1,
                      lambda = c(0, 1, 0),
                      delta = -1)
## criterion value is 1e+24 which will be returned when variance for estimating ED50 is comutationaly negative!
## if we change the tolerance for finding  Moore-Penrose Matrix Inverse to .Machine$double.eps
# when get 2.201179 for the criterion value

equivalence_multiple (x = c(-6.910, -4.6150000, -4.6000000, 6.910),
                      w =   c(0.499998, 0.2230491, 0.2305728, 0.04637965),
                      lx = log(.001), ux = log(1000),
                      param = Theta1,
                      lambda = c(0, 0, 1),
                      delta = -1)

## now let set the real value of the smallest and the largest design point! why?
equivalence_multiple (x = c(log(.001), -4.6150000, -4.6000000, log(1000)),
                      w =   c(0.499998, 0.2230491, 0.2305728, 0.04637965),
                      lx = log(.001), ux = log(1000),
                      param = Theta1,
                      lambda = c(0, 0, 1),
                      delta = -1)

#######################################################################################
######################################################################################

## exponential model. examples for locally optimal design, minimax and maximin
test <- mica(fimfunc = "exp_2par", lx = 0, ux = 1, lp = c(1, 1), up = c(1, 5), iter = 100, k = 3, type = "standardized",
             control = list(seed = 215))

test <- mica(fimfunc = "exp_2par", lx = 0, ux = 1, lp = c(1, 1), up = c(1, 5), iter = 300, k = 2, type = "minimax",
             control = list(seed = 215))
plot(test)



test<- mica(fimfunc = "exp_2par", lx = 0, ux = 1, lp = c(2, 3), up = c(2, 3), iter = 40, k = 2, type = "locally",
            control = list(seed = 215))
plot(test)

mica(fimfunc = "exp_2par", lx = 0, ux = 1, lp = c(2, 3), up = c(2, 3),
     iter = 30, k = 2, type = "locally", control = list(seed = 215, equal_weight = TRUE))


mica(fimfunc = "exp_2par", lx = 0, ux = 1, lp = c(2, 3), up = c(2, 3), iter = 40, k = 2, type = "locally",
     control = list(seed = 215))



test <- mica(fimfunc = "exp_2par", lx = 1, ux = 5, lp = c(2, 3), up = c(2, 3), iter = 500, k = 3, type = "locally",
             control = list(rseed = 215, trace = TRUE, inner_maxiter = 600, ncount = 20, stop_rule = "equivalence", equivalence_every = 100))

test <- mica(fimfunc = "exp_2par", lx = 1, ux = 5, lp = c(2, 3), up = c(2, 3), iter = 500, k = 2, type = "locally",
             control = list(rseed = 215, inner_maxiter = 600,
                            n.seg = 100, ncount = 50,  stop_rule = "equivalence", equal_weight = TRUE))

## equal weight with initial
test <- mica(fimfunc = "exp_2par", lx = 1, ux = 5, lp = c(2, 3), up = c(2, 3), iter = 500, k = 2, type = "locally",
             initial = c(-1, 2, 3, 4),
             control = list(rseed = 215, inner_maxiter = 600,
                            n.seg = 100, ncount = 50,  stop_rule = "equivalence", equal_weight = TRUE))


test <- mica(fimfunc = "exp_2par", lx = 1, ux = 5, lp = c(2, 3), up = c(2, 3), iter = 500, k = 3, type = "locally",
             control = list(rseed = 215, inner_maxiter = 600,
                            n.seg = 100, ncount = 50,  stop_rule = "equivalence", equal_weight = FALSE))

### lgostoc

time1 <- proc.time()
test <- mica(fimfunc = "FIM_logistic",
             lx = -5, ux = 5, lp = c(2, 3), up = c(2, 3), iter = 40, k = 2, type = "locally",
             control = list(seed = 215, trace = T))
time1 <- proc.time() - time1
plot(test)
iterate(test, 10)
#### logistic model examples for symmetric design
time1 <- proc.time()
test <- mica(fimfunc = "logistic", lx = -5, ux = 5, lp = c(0, 1), up = c(3.5 , 1.25), iter = 100, k = 3,
             control = list(rseed = 215, plot_cost = TRUE),type = "minimax")
time1 <- proc.time() - time1

time1 <- proc.time()
test <- mica(fimfunc = "logistic", lx = -5, ux = 5, lp = c(0, 1), up = c(3.5 , 3.5), iter = 100, k = 6,
             control = list(rseed = 215, check_inner_maxeval = FALSE, sym = TRUE, sym_point = (3.5 + 0)/2),
             type = "standardized")
time1 <- proc.time() - time1
plot(test)

test <- mica(fimfunc = "FIM_logistic", lx = -5, ux = 5, lp = c(0, 1), up = c(3.5 , 1.25), iter = 100, k = 5,
             control = list(rseed = 215, sym = TRUE, sym_point = (0 + 3.5)/2),type = "minimax")


test <- mica(fimfunc = "logistic", lx = -5, ux = 5, lp = c(0, 1), up = c(3.5 , 1.25), iter = 100, k = 5,
             control = list(rseed = 215, sym = TRUE, sym_point = (0 + 3.5)/2),type = "minimax")



test <- mica(fimfunc = "logistic", lx = -5, ux = 5,
             lp = c(-1, 1), up = c(1 , 2),
             iter = 100, k = 3,
             type = "minimax",
             control = list(rseed = 215))
plot(test)

#### Enzyme kinetic models. examples for3D plots
mica(fimfunc = "comp_inhibition", lx = c(0, 0), ux = c(30, 60), lp = c(7, 4, 2), up = c(7, 5, 3), k =3, type = "standardized",
     iter = 300, control = list(rseed = 215, inner_maxit = 300, stop_rule = "equivalence",
                                countries = 100, nimperialists = 10))

mica(fimfunc = "comp_inhibition", lx = c(0, 0), ux = c(30, 60), lp = c(7, 4, 2), up = c(7, 5, 3), k =3, type = "standardized",
     iter = 1000, control = list(rseed = 215, inner_space = "vertices", stop_rule = "equivalence",
                                 countries = 100, nimperialists = 10))
## or we can create the set ourselves by lp[1] = up[1], so the fisrt column of param_set was equal to '7'.
param_set <- matrix(c(7, 7, 7, 7, 7, 7, 7, 7, 4, 4, 5, 5, 4, 4, 5, 5, 2, 2, 2, 2, 3, 3, 3, 3), ncol = 3, nrow = 8)
mica(fimfunc = "comp_inhibition", lx = c(0, 0), ux = c(30, 60), lp = c(7, 4, 2), up = c(7, 5, 3), k =3, type = "standardized",
     iter = 1000, control = list(rseed = 215, inner_space = "discrete", stop_rule = "equivalence",
                                 countries = 100, nimperialists = 10, param_set = param_set))


test <- mica(fimfunc = "comp_inhibition", lx = c(0, 0), ux = c(30, 60), lp = c(7, 4, 2), up = c(7, 4, 2), k =3, type = "locally",
             iter = 1000, control = list(rseed = 215,  stop_rule = "equivalence",
                                         countries = 100, equal_weight = TRUE, nimperialists = 10))


#### multi objective optimal design
param <- c(1.563, 1.690, 5.322, 0.137)
## locally optimal design for 4pl model
test <- mica(fimfunc = "FIM_logistic_4par", lx = log(.001) , ux = log(1000),  lp = param,
             up = param, iter = 100, k = 4, type = "locally")

iterate(test, 100)

param <- c(1.563, 1.690, 5.322, 0.137)
time1 <- proc.time()
test <- mica(fimfunc = "logistic_4par", lx = log(.001) , ux = log(1000),  lp = param,
             up = param, iter = 200, k = 4, type = "multiple_locally",
             multiple = list(lambda = rep(1/3, 3), delta = -1))
time1 <- proc.time() - time1




test <- mica(fimfunc = "logistic", lx = -5, ux = 5, lp = c(0, 1), up = c(3.5 , 3.5), iter = 1000, k = 6,
             control = list(rseed = 215, inner_space = "discrete"),
             type = "minimax")


# 2PL model
test <- mica(fimfunc = "FIM_logistic", lx = -3, ux = 3, lp = c(-.1, 2), up = c(1 , 2.5), iter = 100, k = 3,
             control = list(rseed = 215), type = "minimax")





###################### competitive inhibition model
mica (    fimfunc = "FIM_comp_inhibition", lx = c(0, 0),  ux = c(30, 60),
          lp = c(7, 4, 2),
          up = c(7, 5, 3),
          iter = 1000,
          control = list(seed = 1366,
                         inner_maxiter = 300,
                         countries = 100,
                         nimperialists = 10),

          k = 3,
          type = "standardized")




mica (    fimfunc = "FIM_noncomp_inhibition", lx = c(0, 0),  ux = c(30, 60),
          lp = c(7, 4, 2),
          up = c(7, 5, 3),
          iter = 1000,
          control = list(seed = 1366,
                         inner_maxiter = 300,
                         countries = 100,
                         nimperialists = 10),

          k = 3,
          type = "standardized")


mica (    fimfunc = "FIM_uncomp_inhibition", lx = c(0, 0),  ux = c(30, 60),
          lp = c(7, 4, 2),
          up = c(7, 5, 3),
          iter = 1000,
          control = list(seed = 1366,
                         inner_maxiter = 300,
                         countries = 100,
                         nimperialists = 10),

          k = 3,
          type = "standardized")

test <- mica (    fimfunc = "FIM_mixed_inhibition", lx = c(0, 0),  ux = c(30, 60),
                  lp = c(7, 4, 2, 4),
                  up = c(7, 5, 3, 5),
                  iter = 200,
                  initial =   c(3.4614, 0, 4.2801, 3.1426, 30, 0,30, 4.0373, .25, .25, .25, .25),
                  control = list(seed = 1366,
                                 inner_maxiter = 300,
                                 countries = 100,
                                 nimperialists = 10),
                  k = 4,
                  type = "standardized")


## power logitic model
mica (    fimfunc = "FIM_power-logistic",
          lx = -10,
          ux =10,
          lp = c(0, 1),
          up = c(1, 3),
          iter =1000,
          control = list(seed = 1366,
                         inner_maxiter = 200,
                         countries = 200,
                         nimperialists = 20),
          k = 3,
          type = "minimax",
          s = .2)

### log linear



#locally_det<- function(param, lx, ux, npar, fimfunc, control){
locally_det<- function(param,  auxiliary){
  ux <- 1
  lx <- 0
  xstar <- (ux + param[3]) * (lx + param[3]) * (log(ux + param[3]) - log(lx + param[3]))/(ux - lx) - param[3]
  denominator <- dod::det2(FIM_loglin(x = c(lx, xstar, ux) , w = rep(1/3, 3), param = param) , logarithm = FALSE)
  cat("\nThis function is sent to the third optimzation level!!!")
  return(denominator)
}



mica (    fimfunc = "FIM_loglin",
          lx = 0,
          ux =1,
          lp = c(0,1, 1),
          up = c(0,1, 35),
          iter =100,
          control = list(seed = 1366,
                         inner_maxiter = 200,
                         countries = 60,
                         nimperialists = 6),
          k = 4,
          locally = locally_det,
          type = "standardized")

mica (    fimfunc = "FIM_loglin",
          lx = 0,
          ux =1,
          lp = c(0,1, 1),
          up = c(0,1, 35),
          iter =100,
          control = list(seed = 1366,
                         inner_maxiter = 200,
                         countries = 60,
                         nimperialists = 6),
          k = 4,
          type = "standardized")

######### multiple

multica_4pl(   lx = log(.001),
               ux = log(1000),
               param = c(1.563, 1.790, 8.442, 0.137),
               k = 4,
               lambda = rep(1/3, 3),
               delta = -1,
               iter = 150,
               control = list(seed = 1366))



############################################### FW
test <- dod:::mfw(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(1, 1), up = c(1, 5),
                  iter = 300, type = "minimax",
                  control = list(seed = 215))


test <- dod:::mfw(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(1, 1), up = c(1, 1),
                  iter = 300, type = "locally",
                  control = list(seed = 215), initial = c(.3, .6, .2, .8))

test <- dod:::mfw(fimfunc = "FIM_comp_inhibition", lx = c(0, 0),  ux = c(30, 60),
                  lp = c(7, 4, 2),
                  up = c(7, 5, 3),
                  iter = 1000,
                  control = list(seed = 1366),
                  type = "standardized")
test <- iterate(test, 10)


test <- dod:::mfw(fimfunc = "FIM_comp_inhibition", lx = c(0, 0),  ux = c(30, 60),
                  lp = c(7, 4, 2),
                  up = c(7, 4, 2),
                  iter = 50,
                  control = list(seed = 1366),
                  type = "locally")



test_locally <- mica(fimfunc = "FIM_comp_inhibition", lx = c(0, 0),  ux = c(30, 60),
                     lp = c(7, 4, 2),
                     up = c(7, 4, 2),
                     iter = 500,
                     k = 3,
                     type = "locally")


test_locally2 <- mica(fimfunc = arg$FIM,
                      lx =  c(0, 0),
                      ux =  c(30, 60),
                      k = 3,
                      lp = c(7, 4, 2),
                      up = c(7, 4, 2),
                      type = "locally",
                      iter = 500,
                      control = list(trace = TRUE,  plot_cost = FALSE,
                                     stop_rule = "maxiter", plot_deriv = FALSE,
                                     check_every = FALSE, ncount = 100, nimp = 10))
FIM_comp_inhibition_x

test_locally2$evol[[500]]$x

test_locally$evol[[500]]$x





################ for exponential table 4
#exp_res <- vector("list", 4)
i <- 2
time1 <- proc.time()
exp_res[[i]] <- dod:::mfwa(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(1, 1), up = c(1, 5),
                           iter = 2000, type = "standardized",
                           control = list(seed = 215, wtol = .03,n.seg = 10,
                                          stop_rule = "equivalence", stoptol = .998, wtol = .03))
exp_res[[i]]$time1 <- proc.time()-time1

########################################## loglinear

################ forloglinear table 4
#loglin_res <- vector("list", 4)
i <- 3
time1 <- proc.time()
loglin_res[[i]] <- dod:::mfwa(fimfunc = "FIM_loglin", lx = 0, ux = 150, lp = c(1, 1, 1), up = c(1, 1, 35),
                              iter = 200, type = "standardized",
                              control = list(seed = 215, n.seg = 30,
                                             stop_rule = "equivalence", stoptol = .998))
loglin_res[[i]]$time1 <- proc.time()-time1
plot(loglin_res[[i]])


dod:::mfwa(fimfunc = "FIM_loglin", lx = 0, ux = 150, lp = c(1, 1, 1), up = c(1, 1, 35),
           iter = 200, type = "standardized",
           control = list(seed = 215, n.seg = 30,
                          stop_rule = "equivalence", stoptol = .998), initial = )
len <- length(loglin_res[[i]]$best$ncde_x)
ini <- c(loglin_res[[i]]$best$ncde_x, rep(1/len, len))
test <- dod:::mfwa(fimfunc = "FIM_loglin", lx = 0, ux = 150, lp = c(1, 1, 1), up = c(1, 1, 35),
                   iter = 500, type = "standardized",
                   control = list(seed = 215, n.seg = 30,
                                  stop_rule = "equivalence", stoptol = .998), initial = ini)

len <- length(test$best$ncde_x)
ini <- c(test$best$ncde_x, rep(1/len, len))
plot(test)



test <- dod:::mfwa(fimfunc = "FIM_klimpel", lx = 0, ux = 100, lp = c(1, 1), up = c(1, 1),
                   iter = 500, type = "locally",
                   control = list(seed = 215, n.seg = 30,
                                  stop_rule = "equivalence", stoptol = .998))


test <- dod:::mica(fimfunc = "FIM_compartmental_4par", lx = 0, ux = 10, lp = c(1:4), up = c(1:4),
                   k = 4,
                   iter = 100, type = "locally",
                   control = list(seed = 215))


test <- dod:::mica(fimfunc = "FIM_chemical_kinetic", lx = 0.001, ux = 10, lp = c(1:2), up = c(1:2),
                   k = 2,
                   iter = 40, type = "locally",
                   control = list(seed = 215))

test <- dod:::mica(fimfunc = "FIM_klimpel", lx = 0.001, ux = 100, lp = c(1, 1), up = c(1, 1),
                   k = 2,
                   iter = 40, type = "locally",
                   control = list(seed = 215))


det(FIM_compartmental_4par(c(1, 2), c(.5, .5), c(1, 2, 3, 4)))


#################################################################
## H algorithm

dod:::H (fimfunc = "FIM_logistic",
         lx =  -5,
         ux = 5,
         lp = c(.5, 4),
         up = c(1.5, 5),
         iter = 10,
         type = "minimax_D")

dod:::H (fimfunc = "FIM_power_logistic",
         lx =  -5,
         ux = 5,
         lp = c(0, 1),
         up = c(1, 3),
         iter = 10,
         type = "minimax_D", s = .2)


dod:::H (fimfunc = "FIM_power_logistic",
         lx =  -5,
         ux = 5,
         lp = c(.5, 4),
         up = c(1.5, 5),
         iter = 10,
         type = "minimax_D", s = .2)


############################### creating PDF
pack <- "ICAOD"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD", "Rd2pdf", shQuote(path)))

