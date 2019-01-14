#################################
# Exponential growth model
################################
# See how we set stopping rule by adjusting 'stop_rule', 'checkfreq' and 'stoptol'
# It calls the 'senslocally' function every checkfreq = 50 iterations to
# calculate the ELB. if ELB is greater than stoptol = .95, then the algoithm stops.

# initializing by one iteration
res1 <- locally(formula = ~a + exp(-b*x), predvars = "x", parvars = c("a", "b"),
                lx = 0, ux = 1, inipars = c(1, 10),
                iter = 1, k = 2,
                ICA.control= ICA.control(rseed = 100,
                                         stop_rule = "equivalence",
                                         checkfreq = 20, stoptol = .95))
\dontrun{
  # update the algorithm
  res1 <- iterate(res1, 150)
  #stops at iteration 21 because ELB is greater than .95
}

### fixed x, lx and ux are only required for equivalence theorem
\dontrun{
res1.1 <- locally(formula = ~a + exp(-b*x), predvars = "x", parvars = c("a", "b"),
                  lx = 0, ux = 1, inipars = c(1, 10),
                  iter = 100, k = 3,
                  x = c(.25, .5, .75),
                  ICA.control= ICA.control(rseed = 100))
plot(res1.1)
# we can not have an optimal design using this x
}

################################
## two parameter logistic model
################################
res2 <- locally(formula = ~1/(1 + exp(-b *(x - a))),
                predvars = "x", parvars = c("a", "b"),
                family = binomial(), lx = -3, ux = 3,
                inipars = c(1, 3), iter = 1, k = 2,
                ICA.control= list(rseed = 100, stop_rule = "equivalence",
                                  checkfreq = 50, stoptol = .95))
\dontrun{
  res2 <- iterate(res2, 100)
  # stops at iteration 51
}




################################
# A model with two predictors
################################
# mixed inhibition model
\dontrun{
  res3 <- locally(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
                  predvars = c("S", "I"),
                  parvars = c("V", "Km", "Kic", "Kiu"),
                  family = gaussian(),
                  lx = c(0, 0), ux = c(30, 60),
                  k = 4,
                  iter = 300,
                  inipars = c(1.5, 5.2, 3.4, 5.6),
                  ICA.control= list(rseed = 100, stop_rule = "equivalence",
                                    checkfreq = 50, stoptol = .95))
  # stops at iteration 100
}


\dontrun{
  # fixed x
  res3.1 <- locally(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
                  predvars = c("S", "I"),
                  parvars = c("V", "Km", "Kic", "Kiu"),
                  family = gaussian(),
                  lx = c(0, 0), ux = c(30, 60),
                  k = 5,
                  iter = 100,
                  x = c(20, 4, 20, 4, 10,  0, 0, 30, 3, 2),
                  inipars = c(1.5, 5.2, 3.4, 5.6),
                  ICA.control= list(rseed = 100))
}
