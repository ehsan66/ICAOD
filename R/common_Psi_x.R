



Psi_x <- function(x1, mu, FIM,  x, w,  answering){

  ## here x is a degenerate design that putt all its mass on x.
  # x1 is one point
  # This function is required to check the equivalence theorem by ploting and also find the D-efficiency lower bound
  # mu is a vector of measure. For locally optimal design mu = 1
  # FIM: is the Fisher information matrix 
  # x: vector of design points
  # w: vector of design weights
  # answering: the matrix of elements of answering set. Each row is an element.
  # answering set are '\mu'in (3) in ICA paper.
  # return the value of left hand side of (4) in ICA paper as a function of x.
  # if mu = 1  and answering has only one row, then we computing the equivalence theorem left hand side for locally D_optimal design.
  # also can be used for multidimensional models like enzyme kinetic


  if(length(mu) != dim(answering)[1])
    stop("The number of measures is not equal to the number of elements of answering set.")
  if(typeof(FIM) != "closure")
    stop("'FIM' must be of type 'closure'.")



  CardOfRegion <- dim(answering)[2] ## cardinality of region of uncertainty
  n_mu <- dim(answering)[1]


  ## each row is  tr(M^{-1}(\xi, mu_j) %*% I(x, \mu_j))
  # so sum of each row of 'Psi_Point_answering' is \int tr(M^{-1}(\xi, mu) %*% I(x_i, \mu))
  Psi_Point_answering <- matrix(NA, 1,  n_mu)

  i <- 1 ## we need Psi at only point x, so we have only one row
  for(j in 1:n_mu){
    Psi_Point_answering[i,j]  <-
      mu[j] *
      sum(diag(solve(FIM(x = x, w = w, par = answering[j, ]), tol = .Machine$double.xmin) %*%
                 FIM(x = x1, w = 1, par = answering[j, ])
      ))
  }

  PsiAtEachPoint <- rowSums(Psi_Point_answering) - CardOfRegion

  PsiFunction <- PsiAtEachPoint[1]
  return(PsiFunction)
}

