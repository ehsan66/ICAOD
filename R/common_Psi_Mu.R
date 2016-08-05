





Psi_Mu <- function(mu, FIM,  x, w,  answering, PenaltyCoeff){
  # mu is a vector of measure that usually Psi_mu is optimized with respect to.
  # FIM: is the Fisher information matrix as function
  # x: vector of design points
  # w: vector of design weights
  # answering: the matrix of elements of answering set. Each row is an element.
  # answering set are '\mu'in (3) in ICA paper.
  # return the value of left hand side of (4) in ICA paper.
  # if mu = 1  and answering has only one row, then we computing the equivalence theorem left hand side for locally D_optimal design.
 # also can be used for multidimensional models like enzyme kinetic

  if(length(mu) != dim(answering)[1])
    stop("The number of measures is not equal to the number of elements of answering set.")
  if(typeof(FIM) != "closure")
    stop("'FIM' must be of type 'closure'.")



  CardOfRegion <- dim(answering)[2] ##cardinal of region of uncertainty
  n_mu <- dim(answering)[1] ## number of mu, measures


  n_independent <- length(x)/length(w) # number of independent variables.
  one_point_mat <- matrix(x, length(w), n_independent)

  Psi_Point_answering <- matrix(NA, length(w),  n_mu)


  ## each row is  tr(M^{-1}(\xi, mu_j) %*% I(x_i, \mu_j))
  # so sum of each row of 'Psi_Point_answering' is \int tr(M^{-1}(\xi, mu) %*% I(x_i, \mu))
  for(i in 1:length(w)){
    for(j in 1:n_mu){
      Psi_Point_answering[i, j]  <-
        mu[j] *
        sum(diag(solve(FIM(x = x, w = w, par = answering[j, ]), tol = .Machine$double.xmin) %*%
                   FIM(x = as.vector(one_point_mat[i,]), w = 1, par = answering[j, ])
        ))
    }
  }
  #Psi at each points x = poi
  #Psi at each Point
  PsiAtEachPoint <- round(rowSums(Psi_Point_answering) - CardOfRegion, 5)



  ##now we should creat the penalty function. Psi must be zero at each supoport point
  PsiEqualityPenalty <- sum(PenaltyCoeff*pmax(PsiAtEachPoint, 0)^2 + PenaltyCoeff*pmax(-PsiAtEachPoint, 0)^2)
  PsiFunction <-  PsiEqualityPenalty + PenaltyCoeff *  (sum(mu) -1)^2

  return(PsiFunction)
}

