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


## The code for sensitivity analysis in Section 5.
# Results are reported in Table 5

lb <- c(4, 11, 100, 5)
ub <- c(10, 18, 180, 11)

res <- vector("list", 9)
maxEval_vec <- c(1000, 50000)
tol_vec <- c(1e-4, 1e-5)


count <- 1
for (i in 1:length(maxEval_vec)){
  for(j in 1:length(tol_vec)){
    res[[count]] <- bayes(fimfunc = FIM_sig_emax,
                          # formula = ~ theta1 + (theta2 - theta1)*(x^theta4)/(x^theta4 + theta3^theta4),
                          #                 predvars = c("x"), parvars = c("theta1", "theta2", "theta3", "theta4"),
                          lx = .001, ux = 500, k = 7,
                          iter = 500, prior = uniform(lb, ub),
                          ICA.control = list(rseed = 13, ncount = 300, nimp = 30),
                          crt.bayes.control = list(cubature = list(maxEval = maxEval_vec[i],
                                                                   tol = tol_vec[j])))
    count <- count + 1
  }
}


for(l in 1:4){
  tun <- res[[l]]$arg$crt.bayes.control$cubature
  name <- paste(tun$tol, "_", tun$maxEval, ".pdf", sep = "")
  pdf(file = name)
  plot(res[[l]],
       sens.bayes.control = list(
         cubature = list(maxEval = 50000, tol = 1e-6)))
  dev.off()
}



