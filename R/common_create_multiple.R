# model = "logistic_4par"
# type = "multi_locally"
# type = "multi_minimax"

create_multiple <-function(model, fimfunc2, type, ...){
  ## model: a character shows the model
  ## return

  temp <- list(...)
  multi_arg <- temp$multiple

  if (!type %in% c("multiple_locally", "multiple_minimax") )
    stop("'type' must be 'multiple_locally' or 'multiple_minimax'")


  ## we should define this because for very small numbers near zero we have NaN, instead we return -Inf
  log2 <- function(x){
    out <- suppressWarnings(log(x))
    if (is.na(out))
      out <- -1e24
    return(out)
  }




  if (model == "FIM_logistic_4par"){


    if (type == "multiple_locally"){

      lambda <- multi_arg$lambda
      delta <- multi_arg$delta
      s <- 4

      ## definign multi_type
      if (lambda[1] == 1)
        multi_type <- "multiple" else
          if (lambda[2] == 1)
            multi_type <- "ED50" else
              if (lambda[3] == 1)
                multi_type <- "MED" else
                  multi_type <- "multiple"

      multi_type <- "multiple"
      if (multi_type == "multiple"){
        crfunc <- function(param, q, n_independent){
          lq <- length(q)
          n_seg <- lq/(n_independent + 1)
          x_ind <- 1:(n_independent * n_seg)
          w_ind <- (x_ind[length(x_ind)]+1):lq
          x <- q[x_ind]
          w <- q[w_ind]

          #           browser()
          # x <- c(-6.907755, -4.664209, -3.925589, 6.907743)
          # w <- c(0.4806476, 0.4081492, 0.06114303, 0.05006018)
          #################################################################################
          ## compute the information matrix and determine whether it is exactly singular
          FIM_val <- fimfunc2(x=x, w = w, param = param)
          logdet_FIM <- det2(FIM_val, logarithm=TRUE)
          ### we calculate the logarithm to be sure that the information matrix is not singular
          if (length(w) < length(param) || logdet_FIM == 1e+24)
            singular_FIM <- TRUE else
              singular_FIM <- FALSE
          #################################################################################


          #################################################################################
          ## calculating the first part in equation of page 262
          if (lambda[1] != 0){
            if (!singular_FIM)
              part1 <- lambda[1]/s * logdet_FIM else{
                return(1e+24)
              }
            #stop ("singular infomration matrix (logdet is NaN)")
          }else
            part1 <- 0

          #################################################################################

          #################################################################################
          ### computing the generalized inverse of Fisher information matrix
          if (lambda[2] != 0 || lambda[3] != 0){
            #if (singular_FIM)
            FIM_inv <- mpginv(FIM_val, tol = sqrt(.Machine$double.xmin)) #else
            #FIM_inv <- solve(FIM_val, tol =sqrt(.Machine$double.eps))
            #det( FIM_inv)
            #MASS::ginv(FIM_val, tol = sqrt(.Machine$double.eps))
          }
          #################################################################################

          #################################################################################
          ## part2 ED50
          if (lambda[2] == 0)
            part2 <- 0 else{
              ED50_prime <- matrix(c(0, param[3]/param[2]^2, -1/param[2], 0), 1, 4)
              var_ED50 <- ED50_prime %*% FIM_inv %*% t(ED50_prime)
              #see equation 4 of Hyun and Wong (2015)
              ## we dont have negative var!!
              if (var_ED50 > 0){
                # if (lambda[2] == 1)
                #   part2 <- var_ED50 else
                    part2 <- lambda[2] * log(var_ED50)
              } else
                return(1e+24) # because we are finding the minimum

            }
          #################################################################################


          #################################################################################
          ## part3 MED50
          if (lambda[3] == 0)
            part3 <- 0 else {
              ## computing the MED prime
              ########################################### MED prime and var_MED
              if (param[2] > 0)
                MED_prime <- matrix(c(-1/((param[1]+delta)*param[2]),
                                      (param[3] - log(-delta/(param[1]+delta)))/param[2]^2,
                                      -1/param[2], 0),1, 4)
              if (param[2] < 0)
                MED_prime <- matrix(c(1/((param[1]-delta)*param[2]),
                                      (param[3] - log((param[1]-delta)/delta))/param[2]^2,
                                      -1/param[2], 0), 1, 4)
              var_MED <- MED_prime %*% FIM_inv %*% t(MED_prime)
              ###########################################
              ## Eq. 5 Hyun and Wong (2015)
              ## we dont have negative var!!


              if (var_MED > 0){
                # if (lambda[3] == 1)
                #   part3 <- var_MED else
                    part3 <- lambda[3] * log(var_MED)
              } else
                return(1e+24) # because we are finding the minimum

            }
          #################################################################################


          # we find the minimum here
          locally_crfunc <- -part1 + part2 + part3


          # if (locally_crfunc < 0.003)
          #   browser()
          return(locally_crfunc)
        }
      }

    }

    # only for 4PL!
    # equation 6 of Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels
    # Hyun and Wong (2015)
    # directional derivative considerations and show that for a fixed lambda,
    # the sensitivity function for the locally multiple-objective
    d_multi_x <- function(x1, FIM, x, w, param, lambda, delta){

      npar <- 4 #  four parameter logistic

      #################### find the generalized inverse
      FIM_val = FIM(x = x, w = w, param = param)
      #FIM_inv <- solve(FIM_val, tol = .Machine$double.xmin)
      # FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.xmin))
      FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.eps))
      ##################################################


      ####################### writign g
      constant1 <- exp(param[2]*x1 + param[3])
      constant2 <- 1/(1+constant1)
      g_x <- matrix(c(constant2, - param[1] * x1 * constant1*constant2^2, -param[1] * constant1 * constant2^2, 1), 1, 4)
      ##################################################


      ##ED50 prime!
      ED50_prime <- matrix(c(0, param[3]/param[2]^2, -1/param[2], 0), 1, 4)
      ##MED prime!
      if(param[2]>0)
        MED_prime <- matrix(c(-1/((param[1] + delta) * param[2]),
                              (param[3] - log(-delta/(param[1] + delta)))/param[2]^2,
                              -1/param[2], 0), 1, 4)
      if(param[2]<0)
        MED_prime <- matrix(c(1/((param[1]-delta) * param[2]),
                              (param[3] - log((param[1] - delta)/delta))/param[2]^2,
                              -1/param[2], 0), 1, 4)

      ## finally we can write d(x, xi)
      if (lambda[2] == 1){
        d_x_xi <-  g_x %*% FIM_inv %*% t(ED50_prime)  %*% ED50_prime %*% FIM_inv %*% t(g_x) - (ED50_prime %*% FIM_inv %*% t(ED50_prime))
        }else
       # d_x_xi <-  ((g_x %*% FIM_inv %*% t(ED50_prime))^2) - (ED50_prime %*% FIM_inv %*% t(ED50_prime)) else
          if (lambda[3] == 1){
            d_x_xi<-  g_x %*% FIM_inv %*% t(MED_prime)  %*% MED_prime %*% FIM_inv %*% t(g_x) - (MED_prime %*% FIM_inv %*% t(MED_prime))
            #d_x_xi <-  ((g_x %*% FIM_inv %*% t(MED_prime))^2) - (MED_prime %*% FIM_inv %*% t(MED_prime))
      }else{

      d_x_xi <- lambda[1]/npar * g_x %*% FIM_inv %*% t(g_x) +
        lambda[2] * ((g_x %*% FIM_inv %*% t(ED50_prime))^2)/(ED50_prime %*% FIM_inv %*% t(ED50_prime)) +
        lambda[3] * ((g_x %*% FIM_inv %*% t(MED_prime))^2)/(MED_prime %*% FIM_inv %*% t(MED_prime)) - 1
            }
      return(d_x_xi)
    }


  }


  ## lambda and delta are from the global enironmen
  #PsiMulti_x <- function(x1, FIM, x, w, mu, answering, lambda, delta){
  PsiMulti_x <- function(x1, FIM, x, w, mu, answering){

    if(length(mu) != dim(answering)[1])
      stop("The number of measures is not equal to the number of elements of answering set.")

    if(typeof(FIM) != "closure")
      stop("'FIM' must be of type 'closure'.")


    CardOfRegion <- dim(answering)[2] ## cardinality of region of uncertainty
    n_mu <- dim(answering)[1]

    ## each row is  tr(M^{-1}(\xi, mu_j) %*% I(x, \mu_j))
    # so sum of each row of 'Psi_Point_answering' is \int tr(M^{-1}(\xi, mu) %*% I(x_i, \mu))
    Psi_Point_answering <- matrix(NA, 1,  n_mu)


    ##the sum of eah row minus number of parameter is psi
    i <- 1
    for(j in 1:n_mu){
      Psi_Point_answering[i, j]  <-
        mu[j] *
        d_multi_x(x1 = x1, FIM = FIM, x = x, w = w, param = answering[j, ], lambda = lambda, delta = delta)
    }

    PsiAtEachPoint <- rowSums(Psi_Point_answering) #- CardOfRegion
    #x can be any support points, p783 King and Wong(2004). So we choose the first one
    PsiFunction <- PsiAtEachPoint[1]
    return(PsiFunction)

  }

  ## lambda and delta are from the global enironmen
  PsiMulti_Mu <- function(mu, FIM,  x, w, answering, PenaltyCoeff){
    if(length(mu) != dim(answering)[1])
      stop("The number of measures is not equal to the number of elements of answering set.")
    if(typeof(FIM) != "closure")
      stop("'FIM' must be of type 'closure'.")


    CardOfRegion <- dim(answering)[2] ##cardinal of region of uncertainty
    n_mu <- dim(answering)[1] ## number of mu, measures


    n_independent <- length(x)/length(w) # number of independent variables.
    one_point_mat <- matrix(x, length(w), n_independent)

    ##The value of Psi at each design x1 = points[i] and each element of answering
    Psi_Point_answering <- matrix(NA, length(w),  n_mu)

    ##the sum of eah row minus number of parameter is psi
    for(i in 1:length(w)){
      for(j in 1:n_mu){
        Psi_Point_answering[i, j]  <-
          mu[j] *
          d_multi_x(x1 = as.vector(one_point_mat[i,]), FIM = FIM, x = x, w = w,
                    param = answering[j, ], lambda = lambda,  delta = delta)
      }
    }
    #Psi at each points x = poi
    #Psi at each point
    PsiAtEachPoint <- round(rowSums(Psi_Point_answering) - CardOfRegion, 5)
    ##the value of Psi as each design point of design m(x,nu)

    # x can be any support points, p783 King and Wong(2004). So we choose the first one
    PsiEqualityPenalty <- sum(PenaltyCoeff*pmax(PsiAtEachPoint, 0)^2 + PenaltyCoeff*pmax(-PsiAtEachPoint, 0)^2)
    PsiFunction <-  PsiEqualityPenalty + PenaltyCoeff *  (sum(mu) -1)^2

    return(PsiFunction)
  }


  # we dont need psi_Mu for locally version but we kepp it
  return(list(crfunc = crfunc, PsiMulti_Mu = PsiMulti_Mu, PsiMulti_x = PsiMulti_x))
}





