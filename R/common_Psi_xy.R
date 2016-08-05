


## this function is used for plotting the equivalence theorem equaltion for model with two independent variables.
# the function is exactly as psy_x, only with two arument.
# 'Point' will be handeled by the FIM function of the models itself, see 'common_mulit_dimensional_design.R'

Psi_xy <- function(x1, y1, mu, FIM,  x, w,  answering){
  ## WARNINGS: do not change names of 'x1' and 'y1' here unless you check the vectorize in 'PlotPsi_x'
  ## there we have 'Vectorize(FUN = Psi_x, vectorize.args=c("x1", "y1"))'
  ##here x is a degenerate matrix that putt all its mass on x.
  if(length(mu) != dim(answering)[1])
    stop("The number of measures is not equal to the number of elements of answering set.")
  if(typeof(FIM) != "closure")
    stop("'FIM' must be of type 'closure'.")



  CardOfRegion <- dim(answering)[2] # dimension of the region of uncertainty
  n_mu <- dim(answering)[1]

  Psi_Point_answering <- matrix(NA, 1,  n_mu)
  i <- 1
    for(j in 1:n_mu){
      Psi_Point_answering[i,j]  <-
        mu[j] *
        sum(diag(solve(FIM(x = x, w = w, par = answering[j, ]), tol = .Machine$double.xmin) %*%
                           FIM(x = c(x1, y1), w = 1, par = answering[j, ])))
    }

  PsiAtEachPoint <- round(rowSums(Psi_Point_answering) - CardOfRegion, 5)
  PsiFunction <- PsiAtEachPoint[1]

  return(-PsiFunction)
}

