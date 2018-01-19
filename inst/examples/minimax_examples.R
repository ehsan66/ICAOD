########################################
# Two-parameter exponential growth model
########################################
res <- minimax (formula = ~a + exp(-b*x), predvars = "x", parvars = c("a", "b"),
                lx = 0, ux = 1, lp = c(1, 1), up = c(1, 10),
                iter = 1, k = 4, ICA.control= ICA.control(rseed = 100))
\dontrun{
res <- iterate(res, 150)
# iterating the algorithm up to 150 more iterations
}

res # print method
plot(res) # veryfying the general equivalence theorem

########################################
# Two-parameter logistic model.
########################################
# A little bit tickling with the tuning parameters
# reducing the value of maxeval to 200 to increase the speed
cont1 <- crt.minimax.control(optslist = list(maxeval = 500))
cont2 <- ICA.control(rseed = 100, checkfreq = Inf, ncount = 60)

\dontrun{
minimax (formula = ~1/(1 + exp(-b *(x - a))), predvars = "x",
         parvars = c("a", "b"),
         family = binomial(), lx = -3, ux = 3,
         lp = c(0, 1), up = c(1, 2.5), iter = 200, k = 3,
         ICA.control= cont2, crt.minimax.control = cont1)
}

############################################
# An example of a model with two predictors
############################################
# Mixed inhibition model
lower <- c(1, 4, 2, 4)
upper <- c(1, 5, 3, 5)
cont <- crt.minimax.control(optslist = list(maxeval = 100)) # to be faster
res <- minimax(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
               predvars = c("S", "I"),
               parvars = c("V", "Km", "Kic", "Kiu"),
               lx = c(0, 0), ux = c(30, 60), k = 4,
               iter = 1, lp = lower, up = upper,
               ICA.control= list(rseed = 100),
               crt.minimax.control = cont)
\dontrun{
res <- iterate(res, 100)
}
print(res)
plot(res) # sensitivity plot

# Now consider grid points instead of assuming continuous parameter space
# set n.grid to 5
\dontrun{
minimax(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
        predvars = c("S", "I"),
        parvars = c("V", "Km", "Kic", "Kiu"),
        lx = c(0, 0), ux = c(30, 60),
        k = 4, iter = 130, n.grid = 5, lp = lower, up = upper,
        ICA.control= list(rseed = 100, checkfreq = Inf),
        crt.minimax.control = cont)
}

############################################
# Standardized maximin D-optimal designs
############################################
# Now assume the purpose is finding STANDARDIZED designs
# We know from the literature that the locally D-optimal design (LDOD)
# for this model has analytical solution.
# The follwoing function takes the parameter as input and returns
# the design points and weights of LDOD.
# x and w are exactly similar to the arguments of 'fimfunc'.
# x is a vector and returns the design points 'dimension-wise'.
# see explanation of the arguments of 'fimfunc' in 'Details'.

LDOD <- function(V, Km, Kic, Kiu){
  #first dimention is for S and the second one is for I.
  S_min <- 0
  S_max <- 30
  I_min <- 0
  I_max <- 60
  s2 <- max(S_min, S_max*Km*Kiu*(Kic+I_min)/
              (S_max*Kic*I_min+S_max*Kic*Kiu+2*Km*Kiu*I_min+2*Km*Kiu*Kic))
  i3 <- min((2*S_max*Kic*I_min + S_max*Kic*Kiu+2*Km*Kiu*I_min+Km*Kiu*Kic)/
              (Km*Kiu+S_max*Kic), I_max)
  i4 <- min(I_min + (sqrt((Kic+I_min)*(Km*Kic*Kiu+Km*Kiu*I_min+
                                         S_max*Kic*Kiu+S_max*Kic*I_min)/
                            (Km*Kiu+S_max*Kic))), I_max )
  s4 <- max(-Km*Kiu*(Kic+2*I_min-i4)/(Kic*(Kiu+2*I_min-i4)), S_min)
  x <- c(S_max, s2, S_max, s4, I_min, I_min, i3, i4)
  return(list(x = x, w =rep(1/4, 4)))

}
args(LDOD)
\dontrun{
minimax(formula =  ~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu)),
        predvars = c("S", "I"),
        parvars = c("V", "Km", "Kic", "Kiu"),
        lx = c(0, 0), ux = c(30, 60),
        k = 4, iter = 300,
        lp = lower, up = upper,
        ICA.control= list(rseed = 100, checkfreq = Inf),
        crt.minimax.control = cont,
        standardized = TRUE,
        localdes = LDOD)
}


################################################################
# Not necessary!
# The rest of the examples here are only for professional uses.
################################################################
# Imagine you have written your own FIM, say in Rcpp that is faster than
# the FIM created by the formula interface above.

###########################################
# An example of a model with two predictors
###########################################
# For example, th cpp FIM function for the mixed inhibition model is named:
args(FIM_mixed_inhibition)

# We should reparamterize the arguments to match the standard of the
# argument 'fimfunc' (see 'Details').
myfim <- function(x, w, param){
  npoint <- length(x)/2
  S <- x[1:npoint]
  I <- x[(npoint+1):(npoint*2)]
  out <- FIM_mixed_inhibition(S = S, I = I, w = w, param = param)
  return(out)
}
args(myfim)

# Finds minimax optimal design, exactly as before, but NOT using the
# formula interface.
res <- minimax(fimfunc = myfim,
               lx = c(0, 0), ux = c(30, 60), k = 4,
               iter = 1, lp = lower, up = upper,
               ICA.control= list(rseed = 100),
               crt.minimax.control = cont)
\dontrun{
res <- iterate(res, 100)
}
print(res)
plot(res) # sensitivity plot

#########################################
# Standardized maximin D-optimal designs
#########################################
# To match the argument 'localdes' when no formula inteface is used,
# we should reparameterize LDOD.
# The input must be 'param' same as the argument of 'fimfunc'
LDOD2 <- function(param)
  LDOD(V = param[1], Km = param[2], Kic = param[3], Kiu = param[4])

# compare these two:
args(LDOD)
args(LDOD2)
\dontrun{
minimax(fimfunc = myfim,
        lx = c(0, 0), ux = c(30, 60), k = 4,
        iter = 300, lp = lower, up = upper,
        ICA.control= list(rseed = 100, checkfreq = Inf),
        crt.minimax.control = cont,
        standardized = TRUE,
        localdes = LDOD2)
}
