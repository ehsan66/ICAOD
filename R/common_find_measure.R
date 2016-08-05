

find_measure <- function(npar, x, w, answering, FIM, Psi_Mu){
  ## retrun the measure mu for given x, w and A(\xi)
  # FIM is function

if (!is.matrix(answering))
  stop("'answering' must be a matrix")
  if (dim(answering)[1] == 1){
    mu <- 1
    System <- "One equation"
  }else{
    n_mu <-  dim(answering)[1] ## number of  measures

    mu <- optim(
      par = rep(1/dim(answering)[1], dim(answering)[1]),
      fn=Psi_Mu, lower=rep(0, n_mu), upper=rep(1, n_mu),
      control = list(factr =  1e-15),
      method="L-BFGS-B",
      answering = answering,
      x = x,
      w = w,
      FIM = FIM,
      PenaltyCoeff = 5000)$par

#     mu <- directL(fn = Psi_Mu,
#            lower=rep(0, n_mu),
#            upper=rep(1, n_mu),
#            answering = answering,
#            x = x,
#            w = w,
#            FIM = FIM,
#            PenaltyCoeff = 5000,
#             nl.info = FALSE,
#             control=list(xtol_rel=sqrt(.Machine$double.eps),
#                          maxeval=3000))$par
  }

  ##added later
  if (round(sum(mu), 3) != 0)
    mu <- mu/sum(mu)

  return(list(mu = mu))
}

