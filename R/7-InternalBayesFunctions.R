###############################################################################################################*
###############################################################################################################*
plotLLTM_1point <- function(prior, Q, lx, ux, optimalx){


  wlambda <- function(x, w, param, q){
    # x is the design points that are the values for the ability parameters
    #qparam <- apply(q * param, 2, sum)
    qparam <- -sum(q * param)
    # argument of lambda for each param vector
    arg_lambda <- sapply(1:length(qparam), function(k)sum(w * exp(qparam[k] + x)/(1 + exp(qparam[k] + x))^2))
  }

  fim_LLTM <- function(x, w, param){
    Qmat <- Q
    nitems <- nrow(Qmat)
    lmat <- lapply(1:nitems, FUN = function(j)wlambda(q = Qmat[j, ],x = x, w = w, param = param) * Qmat[j, ] %*% t(Qmat[j, ]))
    lmat <- apply(simplify2array(lmat), c(1, 2), sum)
    return(lmat)
  }
  prior_func <- prior$fn
  cr_integrand_D<- function(param, x, w){
    bcrfunc1 <- -1 * det2(fim_LLTM(x = x, w = w, param = param), logarithm= TRUE) * prior_func(t(param))
    dim(bcrfunc1) <- c(1, length(bcrfunc1))
    return(bcrfunc1)
  }
  ability <- seq(from = lx, to = ux, length.out = 100)
  #by = (ux -lx)/50)
  ability <- c(ability, optimalx)

  crtval <- c()
  for (i in 1:length(ability))
    crtval[i] <- hcubature(f = cr_integrand_D, lowerLimit = prior$lower, upperLimit = prior$upper,
                           vectorInterface = FALSE, x = ability[i], w = 1, tol = 1e-3,
                           maxEval = 1000, absError = 0,
                           doChecking = FALSE)$integral

  optimalx_crt <- crtval[length(crtval)]
  crtval <- crtval[order(ability)]
  ability <- ability[order(ability)]

  plot(ability, crtval, type = "l", xlab = expression(theta), ylab = "Bayesian D-criterion")
  points(x = optimalx, y = crtval[length(crtval)],  col = "red" ,pch = 16, cex = 1)
  eff <- optimalx_crt/crtval
  plot(ability, eff, type = "l", xlab = expression(theta), ylab = "Relative Bayesian D-efficiency")
  points(x = optimalx, y =  1,  col = "red" ,pch = 16, cex = 1)
}

# ###############################################################################################################*
# ###############################################################################################################*
control.cubature <- function(tol = 1e-5, maxEval = 50000, absError = 0){
  return(list(tol = tol, maxEval = maxEval, absError = absError))
}
control.quadrature <- function(type = "GLe", level = 8, ndConstruction = "sparse"){
  return(list(type = type, level = level, ndConstruction = ndConstruction))
}
###############################################################################################################*
###############################################################################################################*
Calculate_Cost_bayes <- function(mat, fixed_arg){
  # mat: is the matrix of positions
  # fixed_arg: passed to calculate
  if(!is.matrix(mat))
    stop("'mat' must be matrix")
  n_cost <- dim(mat)[1]
  store <- vector("list", n_cost) #temporarily
  cost <- vector("double", n_cost)
  if(fixed_arg$equal_weight)
    w_equal <- rep(1/fixed_arg$k, fixed_arg$k)
  for(i in 1:dim(mat)[1]){
    x <- mat[i, fixed_arg$x_id]
    if(!fixed_arg$equal_weight)
      w <- mat[i, fixed_arg$w_id] else
        w <- w_equal
      if(fixed_arg$sym){
        x_w <- ICA_extract_x_w(x = x, w = w,
                               sym_point = fixed_arg$sym_point)
        x <- x_w$x
        w <- x_w$w
      }
      store[[i]] <- fixed_arg$crfunc(q = c(x, w), npred = fixed_arg$npred)
  }
  nfeval <-  sum(sapply(store, "[[", "fneval"))
  cost <-  sapply(store, "[[", "val")
  ## we require inner_optima for the ICA functions because they were wriite originally for the minimax designs
  return(list(cost = cost, nfeval = nfeval, inner_optima = matrix(NA, nrow = length(cost))))
}
###############################################################################################################*
###############################################################################################################*

#############################################################################################################*
#############################################################################################################*
create_criterion_bayes <- function(FIM, type, prior, compound, const, multiple, localdes = NULL, method, npar, crt.bayes.control, IRTpars){
  prior_func <- prior$fn
  if (is.null(npar))
    npar <- prior$npar
  lp <- prior$lower
  up <- prior$upper

  if (method == "cubature"){
    #### make the cr_integrand

    if (type == "DPA")
      cr_integrand_DPA <- function(param, x, w){

        if (compound$alpha != 0){
          bcrfunc1 <- -(compound$alpha/4  * sapply(FIM(x = x, w = w, param = t(param)), det2, logarithm= TRUE) +
                          apply(param, 2, function(col_par)(1 -compound$alpha) * log(sum(w * compound$prob(x = x, param = col_par))))) * prior_func(t(param))
        }else{
          bcrfunc1 <- -(apply(param, 2, function(col_par)(1 -compound$alpha) * log(sum(w * compound$prob(x = x, param = col_par))))) * prior_func(t(param))

        }
        # if (any(is.infinite(bcrfunc1)))
        #   bcrfunc1[(is.infinite(bcrfunc1))] <- 0
        dim(bcrfunc1) <- c(1, length(bcrfunc1))
        return(bcrfunc1)
      }
    if (type == "DPM")

      cr_integrand_DPM <- function(param, x, w){
        bcrfunc1 <- apply(param, 2,  FUN = function(col_par)-(compound$alpha/4 * det2(FIM(x = x, w = w, param = col_par), logarithm = TRUE) + (1 -compound$alpha) * log(min(compound$prob(x = x, param = col_par))))) * prior_func(t(param))
        if(any(is.infinite(bcrfunc1)))
          bcrfunc1[(is.infinite(bcrfunc1))] <- 0
        dim(bcrfunc1) <- c(1, length(bcrfunc1))
        return(bcrfunc1)
      }
    if (type == "D")
      cr_integrand_D<- function(param, x, w){
        bcrfunc1 <- -1 * sapply(FIM(x = x, w = w, param = t(param)), det2, logarithm= TRUE) * prior_func(t(param))
        # FIM(x = x, w = w, param = t(param))
        # FIM_logistic(x =x, w = w, param = c(0, 1.05))
        #bcrfunc1 <- apply(param, 2, FUN = function(col_par)-det2(FIM(x = x, w = w, param = col_par), logarithm = TRUE)) * prior_func(t(param))
        dim(bcrfunc1) <- c(1, length(bcrfunc1))
        return(bcrfunc1)
      }

    if (type == "D_LLTM"){
      wlambda <- function(x, w, param, q){
        # x is the design points that are the values for the ability parameters
        qparam <- apply(q * param, 2, sum)
        # argument of lambda for each param vector
        arg_lambda <- sapply(1:length(qparam), function(k)sum(w * exp(qparam[k] + x)/(1 + exp(qparam[k] + x))^2))
      }

      crt_LLTM <- function(x, w, param){
        Qmat <- cbind(1, IRTpars$Q)
        nitems <- nrow(Qmat)
        lmat <- sapply(1:nitems, FUN = function(j) wlambda(q = Qmat[j, ],x = x, w = w, param = param))
        ncolQ <- ncol(Qmat)

        lmat <- sapply(1:nitems, FUN = function(j) ((wlambda(q = Qmat[j, ],x = x, w = w, param = param))^ncolQ) * det(Qmat[j, ] %*% t(Qmat[j, ])))
        #lmat ^ ncol(Qmat) * det(Qmat %*% t(Qmat))


        #lmat <- log(lmat)
        #crtval <- matrix(apply(lmat, 1, sum))
        return(lmat)
      }
    }

    if(type == "multiple"){
      stop("BUG: No Bayesian multiple objective optimal designs is implemented yet!")
      # temp_multiple <- create_multiple (model =  "FIM_logistic_4par", fimfunc2 = FIM, multiple_arg = multiple)
      # crfunc_multiple <- function(param, x, w){
      #   temp_multiple$crfunc(param =param, x =x, w=w, multiple = multiple)
      # }
      # cr_integrand_multiple <- function(param, x, w){
      #   bcrfunc1 <- apply(param, 2, FUN = function(col_par)crfunc_multiple(x = x, w = w, param = col_par)) * prior_func(t(param))
      #   dim(bcrfunc1) <- c(1, length(bcrfunc1))
      #   return(bcrfunc1)
      # }
    }
    cr_integrand <- switch(type, "D" = cr_integrand_D, "DPA" = cr_integrand_DPA, "DPM" = cr_integrand_DPM, "D_LLTM" = crt_LLTM)
    #cr_integrand <- switch(type, "D" = cr_integrand_D, "DPA" = cr_integrand_DPA, "DPM" = cr_integrand_DPM, "multiple" = cr_integrand_multiple)
    crfunc_bayesian  <- function(q, npred) {
      lq <- length(q) # q is the design points and weights
      pieces <- lq / (npred + 1)
      x_ind <- 1:(npred * pieces)
      w_ind <- (x_ind[length(x_ind)] + 1):lq
      x <- q[x_ind]
      w <- q[w_ind]

      out <- hcubature(f = cr_integrand, lowerLimit = lp, upperLimit = up, vectorInterface = TRUE, x = x, w = w,
                       tol = crt.bayes.control$cubature$tol, maxEval = crt.bayes.control$cubature$maxEval,
                       absError = 0,# norm = crt.bayes.control$cubature$norm,
                       doChecking = FALSE)

      if (out$integral == Inf || is.nan(out$integral) || out$integral == -Inf )
        out$integral <- 1e+20
      if (const$use)
        val <- out$integral + 5000 * (sum(w) - 1)^2 + const$pen_func(x) else
          val <- out$integral + 5000 * (sum(w) - 1)^2
      return(list(val = val, fneval = out$functionEvaluations))
    }
    crfunc <- crfunc_bayesian
  }#cubature
  #############################################################################*

  #############################################################################*
  ###quadrature methods
  if (method == "quadrature"){
    # if (type == "D"){
    #   cr_integrand_D<- function(param, x, w){
    #     bcrfunc1 <- -1 * sapply(FIM(x = x, w = w, param = param), det2, logarithm= TRUE) * prior_func(param)
    #     if (dim(param)[1] != dim(bcrfunc1)[1])
    #       bcrfunc1 <- t(bcrfunc1)
    #     return(bcrfunc1)
    #   }
    # }
    #
    #
    # if (type == "DPA")
    #   cr_integrand_DPA <- function(param, x, w){
    #     bcrfunc1 <- -(compound$alpha/4  * sapply(FIM(x = x, w = w, param = param), det2, logarithm= TRUE) +
    #                     apply(param, 1, function(row_par)(1 -compound$alpha) * log(sum(w * compound$prob(x = x, param = row_par))))) * prior_func(param)
    #     if(any(is.infinite(bcrfunc1)))
    #       bcrfunc1[(is.infinite(bcrfunc1))] <- 0
    #     if (dim(param)[1] != dim(bcrfunc1)[1])
    #       bcrfunc1 <- t(bcrfunc1)
    #     return(bcrfunc1)
    #   }
    #
    # cr_integrand <- switch(type, "D" = cr_integrand_D, "DPA" = cr_integrand_DPA, "DPM" = cr_integrand_DPM, "D_LLTM" = crt_LLTM)
    #
    #
    # nw <- createNIGrid(dim = npar, type = crt.bayes.control$quadrature$type,
    #                    level = crt.bayes.control$quadrature$level,
    #                    ndConstruction = crt.bayes.control$quadrature$ndConstruction)
    # rescale(nw, domain = matrix(c(prior$lower, prior$upper), ncol=2))
    #
    # crfunc_bayesian  <- function(q, npred) {
    #   lq <- length(q) # q is the design points and weights
    #   pieces <- lq / (npred + 1)
    #   x_ind <- 1:(npred * pieces)
    #   w_ind <- (x_ind[length(x_ind)] + 1):lq
    #   x <- q[x_ind]
    #   w <- q[w_ind]
    #
    #   out <- quadrature(f = cr_integrand, grid = nw, x=x, w=w)
    #
    #   if (out == Inf || is.nan(out) || out == -Inf )
    #     out <- 1e+20
    #   if (const$use)
    #     val <- out + 5000 * (sum(w) - 1)^2 + const$pen_func(x) else
    #       val <- out + 5000 * (sum(w) - 1)^2
    #   return(list(val = val, fneval = dim(nw$nodes)[1]))
    # }
    # crfunc <- crfunc_bayesian
    stop("quadrature methods are not implemented yet!")
  }




  #####################################################################*
  return(list(crfunc = crfunc))
}
###############################################################################################################*
###############################################################################################################*
create_Psi_bayes <- function(type, prior, FIM, lp, up, npar, truncated_standard, const, sens.bayes.control, compound, IRTpars, method){

  if (method == "cubature"){
    if (type == "D"){
      Psi_x_integrand_D  <- function(x1, param,  FIM, x, w, prior_func){

        # NOTE: in cubature when integrand is vectorized:
        #  1- the input of the integrand (param) is a MATRIX with  length(mean) or length(lp) number of rows and l number of columns  (l is for vectorization)
        #  2- the integrand should return a matrix of 1 * l (vectorization).
        ## param: as matrix: EACH COLUMN is one set of parameters
        ## return a matrix with 1 * length(param) dimension
        FIM_x <- FIM(x = x, w = w, par = t(param))
        FIM_x1 <- FIM(x = x1, w = 1, par = t(param))
        deriv_integrand <- sapply(1:ncol(param), FUN = function(j)sum(diag(solve(FIM_x[[j]])%*%FIM_x1[[j]]))) * prior_func(t(param))

        # deriv_integrand <- apply(param, 2, FUN = function(col_par)sum(diag(solve(FIM(x = x, w = w, par = col_par)) %*%
        #                                                                       FIM(x = x1, w = 1, par = col_par)))) * prior_func(t(param))
        if (is.null(dim(deriv_integrand)))
          dim(deriv_integrand) <- c(1, length(deriv_integrand))
        return(deriv_integrand)
      }

      Psi_x_bayes <- function(x1,  x, w){
        ## here x is a degenerate design that putt all its mass on x.
        # x1 is one point
        # This function is required to check the equivalence theorem by ploting and calculating the D-efficiency lower bound
        # prior funtion
        # FIM: is the Fisher information matrix
        # x: vector of design points
        # w: vector of design weights
        # we have npar here beacuse in optimization we shorten lp and up when the lower bound and uper bound is the same
        # truncated_standard: required for the prior because the prior is truncated
        out <- hcubature(f = Psi_x_integrand_D, lowerLimit = lp, upperLimit = up, vectorInterface = TRUE,
                         x = x, w = w, prior_func = prior$fn, FIM = FIM, x1 = x1,
                         tol = sens.bayes.control$cubature$tol,
                         maxEval = sens.bayes.control$cubature$maxEval,
                         doChecking = FALSE,
                         #norm = sens.bayes.control$cubature$norm,
                         absError = 0)
        return(out$integral/truncated_standard  - npar)
      }

      Psi_xy_bayes <- function(x1, y1, x, w){


        ## this function is used for plotting the equivalence theorem equation for model with two independent variables.
        # the function is exactly as psy_x_bayes, only with two arument.
        # 'Point' will be handeled by the FIM function of the models itself, see 'common_mulit_dimensional_design.R'
        ## WARNINGS: do not change names of 'x1' and 'y1' here unless you check the vectorize in 'PlotPsi_x'
        ## there we have 'Vectorize(FUN = Psi_x, vectorize.args=c("x1", "y1"))'

        ## here x is a degenerate design that putt all its mass on x.
        # x1 is one point
        # This function is required to check the equivalence theorem by ploting and also find the D-efficiency lower bound
        # prior funtion
        # FIM: is the Fisher information matrix
        # x: vector of design points
        # w: vector of design weights
        # we have npar here beacuse in optimization we shorten lp and up when the lower bound and uper bound is the same
        # truncated_standard: required for the prior because the prior is truncated

        out <- hcubature(f = Psi_x_integrand_D, lowerLimit = lp, upperLimit = up, vectorInterface = TRUE,
                         x = x, w = w, prior_func = prior$fn, FIM = FIM, x1 = c(x1, y1),
                         tol = sens.bayes.control$cubature$tol,
                         maxEval = sens.bayes.control$cubature$maxEval,
                         doChecking = FALSE,
                         #norm = sens.bayes.control$cubature$norm,
                         absError = 0)
        #return(list(val = out$integral/truncated_standard  - npar, nfeval = out$functionEvaluations))
        return(-(out$integral/truncated_standard  - npar))
      }
    }
    if (type == "DPA" || type == "DPAM"){

      Psi_x_integrand_DP <- function(x1,  FIM,  x, w, param, npar, alpha,  type, prob, prior_func){
        # NOTE: in cubature when integrand is vectorized:
        #  1- the input of the integrand (param) is a MATRIX with length(mean) or length(lp) number of rows and l number of columns  (l is for vectorization)
        #  2- the integrand should return a matrix of 1 * l (vectorization).
        ## param: as matrix: EACH COLUMN is one set of parameters
        if (type == "DPA"){

          FIM_x <- FIM(x = x, w = w, par = t(param))
          FIM_x1 <- FIM(x = x1, w = 1, par = t(param))
          if (compound$alpha != 0){
            deriv_integrand <-  sapply(1:ncol(param), FUN = function(j)alpha/npar * sum(diag(solve(FIM_x[[j]]) %*% FIM_x1[[j]])) +
                                         (1-alpha) * (prob(x1, param[, j])- sum(w * prob(x, param[, j])))/sum(w * prob(x, param[, j]))) * prior_func(t(param))
          }else{
            deriv_integrand <-  sapply(1:ncol(param), FUN = function(j) (prob(x1, param[, j])- sum(w * prob(x, param[, j])))/sum(w * prob(x, param[, j]))) * prior_func(t(param))
          }

          # deriv_integrand <- apply(param, 2,
          #                          FUN = function(col_par)alpha/npar * sum(diag(solve(FIM(x = x, w = w, par = col_par)) %*% FIM(x = x1, w = 1, par = col_par)))
          #                          + (1-alpha) * (prob(x1, col_par)- sum(w * prob(x, col_par)))/sum(w * prob(x, col_par))) * prior_func(t(param))
        }else
          if (type == "DPM"){
            deriv_integrand <- apply(param, 2, FUN = function(col_par)alpha/npar * sum(diag(solve(FIM(x = x, w = w, par = col_par)) %*% FIM(x = x1, w = 1, par = col_par))) + (1-alpha) * (prob(x1, col_par)- min(prob(x, col_par)))/min(prob(x, col_par))) * prior_func(t(param))
          }else
            stop("BUG: check 'type'")
        dim(deriv_integrand) <- c(1, length(deriv_integrand))
        return(deriv_integrand)
      }
      Psi_x_bayes <- function(x1, x, w){
        ## here x is a degenerate design that putt all its mass on x.
        # x1 is one point
        # This function is required to check the equivalence theorem by ploting and also find the D-efficiency lower bound
        # prior funtion
        # FIM: is the Fisher information matrix
        # x: vector of design points
        # w: vector of design weights
        # we have npar here beacuse in optimization we shorten lp and up when the lower bound and uper bound is the same
        # truncated_standard: required for the prior because the prior is truncated
        out <- hcubature(f = Psi_x_integrand_DP, lowerLimit = lp, upperLimit = up, vectorInterface = TRUE,
                         x = x, w = w, prior_func = prior$fn, FIM = FIM, x1 = x1,
                         alpha = compound$alpha,  type = type, prob = compound$prob, npar = npar,
                         tol = sens.bayes.control$cubature$tol,
                         maxEval = sens.bayes.control$cubature$maxEval,
                         doChecking = FALSE,
                         #norm = sens.bayes.control$cubature$norm,
                         absError = 0)
        return(out$integral/truncated_standard  - compound$alpha)
      }


      Psi_xy_bayes <- function(x1, y1, x, w){
        ## WARNINGS: do not change names of 'x1' and 'y1' here unless you check the vectorize in 'PlotPsi_x'
        ## there we have 'Vectorize(FUN = Psi_x, vectorize.args=c("x1", "y1"))'
        ## this function is used for plotting the equivalence theorem equation for models with two independent variables.
        # the function is exactly as psy_x_bayes_compound, only with two aruments.
        # 'Point' will be handeled by the FIM function of the models itself, see 'common_mulit_dimensional_design.R'

        # here x is a degenerate design that putt all its mass on x.
        # x1 is one point
        # This function is required to check the equivalence theorem by ploting and also find the D-efficiency lower bound
        # prior funtion
        # FIM: is the Fisher information matrix
        # x: vector of design points
        # w: vector of design weights
        # we have npar here beacuse in optimization we shorten lp and up when the lower bound and uper bound is the same
        # truncated_standard: required for the prior because the prior is truncated

        out <- hcubature(f = Psi_x_integrand_DP, lowerLimit = lp, upperLimit = up, vectorInterface = TRUE,
                         alpha = compound$alpha,  type = type, prob = compound$prob, npar = npar,
                         x = x, w = w, prior_func = prior$fn, FIM = FIM, x1 = c(x1, y1),
                         tol = sens.bayes.control$cubature$tol,
                         maxEval = sens.bayes.control$cubature$maxEval,
                         doChecking = FALSE,
                         #norm = sens.bayes.control$cubature$norm,
                         absError = 0)
        return(-(out$integral/truncated_standard  - compound$alpha))
      }
    }
    if (type == "D_LLTM"){
      M_i <- function(x, w, param, q){
        # x is the design points that are the ability parameter
        temp <-  exp(sum(q * param) + x)
        sum(w * temp/(1 + temp)^2)^(length(q) + 1) * q %*% t(q)
      }


      Psi_x_integrand_D_LLTM  <- function(x1, param, x, w,Q, prior_func){

        browser()
        Q <- cbind(1, IRTpars$Q)
        nitems <- nrow(Q)
        M <- sapply(1:nitems, FUN = function(j) M_i(q = Q[j, ],x = x, w = w, param = param))
        M_x <- sapply(1:nitems, FUN = function(j) M_i(q = Q[j, ],x = x1, w = 1, param = param))

        deriv_integrand <- sapply(1:ncol(param), FUN = function(j)sum(diag(solve(M[[j]])%*%M_x[[j]]))) * prior_func(t(param))
        #deriv_integrand <- sapply(1:ncol(param), FUN = function(j)sum(diag(solve(FIM_x[[j]])%*%FIM_x1[[j]]))) * prior_func(t(param))

        # deriv_integrand <- apply(param, 2, FUN = function(col_par)sum(diag(solve(FIM(x = x, w = w, par = col_par)) %*%
        #                                                                       FIM(x = x1, w = 1, par = col_par)))) * prior_func(t(param))
        if (is.null(dim(deriv_integrand)))
          dim(deriv_integrand) <- c(1, length(deriv_integrand))
        return(deriv_integrand)
      }

      Psi_x_bayes <- function(x1,  x, w){
        ## here x is a degenerate design that putt all its mass on x.
        # x1 is one point
        # This function is required to check the equivalence theorem by ploting and calculating the D-efficiency lower bound
        # prior funtion
        # FIM: is the Fisher information matrix
        # x: vector of design points
        # w: vector of design weights
        # we have npar here beacuse in optimization we shorten lp and up when the lower bound and uper bound is the same
        # truncated_standard: required for the prior because the prior is truncated
        out <- hcubature(f = Psi_x_integrand_D_LLTM, lowerLimit = lp, upperLimit = up, vectorInterface = TRUE,
                         x = x, w = w, prior_func = prior$fn, Q = IRTpars$Q, x1 = x1,
                         tol = sens.bayes.control$cubature$tol,
                         maxEval = sens.bayes.control$cubature$maxEval,
                         doChecking = FALSE,
                         absError = 0)
        return(out$integral/truncated_standard  - npar)
      }
      Psi_xy_bayes <- NA

    }
  }
  if (method == "quadrature"){

    # nw <- createNIGrid(dim = npar, type = sens.bayes.control$quadrature$type,
    #                    level = sens.bayes.control$quadrature$level,
    #                    ndConstruction = sens.bayes.control$quadrature$ndConstruction)
    # rescale(nw, domain = matrix(c(prior$lower, prior$upper), ncol=2))
    #
    # if (type == "D"){
    #   Psi_x_integrand_D  <- function(x1, param,  FIM, x, w, prior_func){
    #
    #     ## param: as matrix: EACH row is one set of parameters
    #     ## return a matrix with 1 * length(param) dimension
    #     FIM_x <- FIM(x = x, w = w, par = param)
    #     FIM_x1 <- FIM(x = x1, w = 1, par = param)
    #     deriv_integrand <- sapply(1:nrow(param), FUN = function(j)sum(diag(solve(FIM_x[[j]])%*%FIM_x1[[j]]))) * prior_func(param)
    #
    #     if (dim(param)[1] != dim(deriv_integrand)[1])
    #       deriv_integrand  <- t(deriv_integrand)
    #     return(deriv_integrand)
    #   }
    #
    #   Psi_x_bayes <- function(x1,  x, w){
    #     ## here x is a degenerate design that putt all its mass on x.
    #     # x1 is one point
    #     # This function is required to check the equivalence theorem by ploting and calculating the D-efficiency lower bound
    #     # prior funtion
    #     # FIM: is the Fisher information matrix
    #     # x: vector of design points
    #     # w: vector of design weights
    #     # we have npar here beacuse in optimization we shorten lp and up when the lower bound and uper bound is the same
    #     # truncated_standard: required for the prior because the prior is truncated
    #
    #     out <- quadrature(f = Psi_x_integrand_D, grid = nw, x=x, w=w, FIM = FIM, x1 = x1, prior_func = prior$fn)
    #     return(out/truncated_standard  - npar)
    #   }
    #
    #   Psi_xy_bayes <- function(x1, y1, x, w){
    #
    #
    #     ## this function is used for plotting the equivalence theorem equation for model with two independent variables.
    #     # the function is exactly as psy_x_bayes, only with two arument.
    #     # 'Point' will be handeled by the FIM function of the models itself, see 'common_mulit_dimensional_design.R'
    #     ## WARNINGS: do not change names of 'x1' and 'y1' here unless you check the vectorize in 'PlotPsi_x'
    #     ## there we have 'Vectorize(FUN = Psi_x, vectorize.args=c("x1", "y1"))'
    #
    #     ## here x is a degenerate design that putt all its mass on x.
    #     # x1 is one point
    #     # This function is required to check the equivalence theorem by ploting and also find the D-efficiency lower bound
    #     # prior funtion
    #     # FIM: is the Fisher information matrix
    #     # x: vector of design points
    #     # w: vector of design weights
    #     # we have npar here beacuse in optimization we shorten lp and up when the lower bound and uper bound is the same
    #     # truncated_standard: required for the prior because the prior is truncated
    #     browser()
    #     out <- quadrature(f = Psi_x_integrand_D, grid = nw, x=x, w=w, FIM = FIM,  x1 = c(x1, y1), prior_func = prior$fn)
    #     return(-(out/truncated_standard  - npar))
    #   }
    # }

  stop("not implemented yet!")

  }
  return(list(Psi_xy_bayes = Psi_xy_bayes, Psi_x_bayes = Psi_x_bayes))
}
######################################################################################################*
######################################################################################################*
PlotPsi_x_bayes <- function(x, w, lower, upper, Psi_x, const, plot_3d  = "lattice"){
  # plot the equivalanece theorem, for both one and two dimensional.
  # x: vector of the point of the design
  # w: vector of the weights of the design
  # lower: lower bound
  # psi_x psi_x as a function of x. in multiobjective optimal design ('multi_ICA'), 'PsiMulti_x' is 'given as Psi_x'
  # FIM: fisher information matrix function
  # prior: the prior function that is necessary for bayesian optimal design
  #: lp and up: the upper bound and lower bound of the parameter space required for the integration
  # vector of the measure. For locally optimal design mu = 1.
  # plot_3d: which package should be used to plot the 3d plots, 'rgl' or 'lattice'
  if(length(lower) == 1){
    xPlot <- sort(c(seq(lower,upper, length.out = ceiling(100 + upper - lower)), x))
    #xPlot <- seq(lower,upper, length.out = 250)
    PsiPlot<- sapply(1:length(xPlot), FUN = function(j) Psi_x(x1 = xPlot[j], x = x, w = w))
    plot(xPlot, PsiPlot, type = "l",
         col = "blue",   xlab = "Design interval",
         ylab = "c", xaxt = "n")
    abline(h = 0, v = c(x) ,col = "grey", lty = 3)

    Point_y <- sapply(1:length(x), FUN = function(j) Psi_x(x1 = x[j], x = x, w = w))

    points(x = x,  y = Point_y, col = "red" ,pch = 16, cex = 1)

    axis(side = 1, at = c(lower, x, upper, 0),
         labels = round(c(lower, x, upper, 0), 4))
  }

  if(length(lower) == 2){
    xlab_3d <- "x1"
    ylab_3d <- "x2"
    len <- 20
    xPlot <- seq(lower[1],upper[1], length.out = len)
    yPlot <- seq(lower[2],upper[2], length.out = len)

    ncolor = 100
    color1 <- rev(rainbow(ncolor, start = 0/6, end = 4/6))

    if (plot_3d  == "lattice"){
      xy <- expand.grid(xPlot, yPlot)
      if (const$use){
        pen_val <- apply(xy, 1, const$pen_func_1point)
        temp_ind <- which(pen_val == 0)
        xy <- xy[temp_ind, ]
        pp <- cbind(x[1:(length(x)/2)], x[(length(x)/2 + 1):length(x)])
        xy <- rbind(as.matrix(xy), pp)
      }
      colnames(xy) <- c("x", "y") ## require for wireframe
      z <- - apply(xy, 1, FUN = function(row1)Psi_x(x1 = row1[1], y1= row1[2], x = x, w = w))
      if (requireNamespace("lattice", quietly = TRUE)){
        wireframe_dat <- as.data.frame(xy)
        wireframe_dat$z <- as.vector(z)
        p1 <- lattice::wireframe( z ~ x * y, data = wireframe_dat,
                                  col.regions=  color1,
                                  drap = TRUE,
                                  scales = list(arrows = FALSE),
                                  shade = FALSE,
                                  xlab = xlab_3d, ylab = ylab_3d,
                                  pretty = TRUE,
                                  zlab = "c",
                                  colorkey = list(col = color1, tick.number = 16))
        print(p1)
        ###contour plot
        xyz <- as.data.frame(xy)
        xyz$z <- as.vector(z)
        p2 <- lattice::contourplot( z ~ x * y, data = xyz,
                                    xlab = xlab_3d,
                                    ylab = ylab_3d,
                                    region = TRUE,
                                    col.regions = color1,
                                    colorkey = list(col = color1, tick.number = 16),
                                    cuts = 13)
        print(p2)
      }else {
        warning("Package 'lattice' is not installed in your system. The 3D derivation plot can not be plotted unless this packages is installed.")
      }
    }
    if(plot_3d  == "rgl"){

      if (requireNamespace("rgl", quietly = TRUE)){
        rgl::.check3d()
        gg  <- function(xx, yy)
          Psi_x(x1 = xx, y1= yy, x = x, w = w)

        z <- outer(xPlot, yPlot, Vectorize(gg))

        if (const$use){
          ff <- function(xx, yy) res <- const$pen_func_1point(c(xx, yy))
          pen_val <- outer(xPlot, yPlot, Vectorize(ff))
          temp_ind <- which(pen_val == 0)
          z[-temp_ind] <- NA
        }
        zcol  <- cut(as.vector(z), ncolor)
        rgl::persp3d(x = xPlot, y = yPlot,
                     z = z,
                     col = color1[zcol],
                     smooth = FALSE,
                     alpha = .8,
                     shininess    = 128,
                     xlab = xlab_3d, ylab = ylab_3d, zlab = "c(x, \\xi, \\mu)")
        Point_mat <- matrix(round(x, 3), ncol = length(lower), nrow = length(x)/length(lower))
        ##adding point
        Point_y <- sapply(1:dim(Point_mat)[1], FUN = function(j) Psi_x(x1 = Point_mat[j, 1],
                                                                       y1= Point_mat[j, 2],
                                                                       x=x,
                                                                       w=w))
        rgl::points3d(x = Point_mat[, 1],
                      y = Point_mat[, 2],
                      z = Point_y,
                      col = "darkred",
                      size = 11,
                      point_antialias = TRUE)
        text_point <- paste("(", xlab_3d, "=", round(Point_mat[, 1], 3), ",", ylab_3d, "=", round(Point_mat[, 2], 3), ";c=",  round(Point_y,3), ")", sep ="")
        rgl::text3d(x = Point_mat[, 1],
                    y = Point_mat[, 2],
                    z = Point_y + .1,
                    texts = text_point,
                    col = "darkred",
                    font = 4)

        rgl::grid3d("z")
      }else{
        warning("Package 'rgl' is not installed in your system. The 3D derivation plot can not be plotted unless this packages is installed.")
      }
    }
  } ## length = 2
}
######################################################################################################*
######################################################################################################*


Psi_x_minus_bayes <- function(x1, x, w, Psi_x_bayes, const){
  ## mu and answering are only to avoid having another if when we want to check the maximum of sensitivity function
  if (const$use)
    Out <- -Psi_x_bayes(x1 = x1, x = x, w = w) + const$pen_func_1point(x1) else
      Out <- -Psi_x_bayes(x1 = x1, x = x, w = w)
    return(Out)
}
######################################################################################################*
######################################################################################################*
bayes_inner <- function(fimfunc = NULL,
                        formula,
                        predvars,
                        parvars,
                        family = gaussian(),
                        lx,
                        ux,
                        type = c("D", "DPA", "DPM", "multiple"),
                        crt_method = c("cubature", "quadrature"),
                        iter,
                        k,
                        npar = NULL,
                        prior = list(),
                        sens_method = c("cubature", "quadrature"),
                        compound = list(prob = NULL, alpha = NULL),
                        multiple.control = list(),
                        ICA.control =  list(),
                        crt.bayes.control = list(),
                        sens.bayes.control = list(),
                        initial = NULL,
                        const = list(ui = NULL, ci = NULL, coef = NULL),
                        plot_3d = c("lattice", "rgl"),
                        IRTpars = NULL,
                        ...) {
  time1 <- proc.time()
  fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =lx, ux = ux,
                                       iter = iter, k = k,
                                       paramvectorized = TRUE, prior = prior)
  if(missing(formula)){
    # to handle ...
    fimfunc2 <- function(x, w, param)
      lapply(1:dim(param)[1], FUN = function(j)fimfunc(x = x, w = w, param = param[j, ]), ...)
    #fimfunc(x = x, w = w, param = param,...)
    fimfunc_sens <- function(x, w, param)
      fimfunc(x = x, w = w, param = param,...) ## it can not be equal to fimfunc2 when inner_space is discrete
  } else{

    #if (length(prior$lower) != length(parvars))
    if (length(prior$lower) != fimfunc_formula$num_unknown_param)
      stop("length of 'prior$lower' is not equal to the number of unknown (not fixed) parameters")
    # fim_localdes <- fimfunc_formula$fimfunc_formula
    fimfunc2 <- fimfunc_formula$fimfunc_formula ## can be vectorized with respect to parameters!
    fimfunc_sens <- fimfunc_formula$fimfunc_sens_formula
  }

  #######################################################*
  ICA.control <- do.call("ICA.control", ICA.control)
  ICA.control <- add_fixed_ICA.control(ICA.control.list = ICA.control)

  sens.bayes.control <- do.call("sens.bayes.control", sens.bayes.control)
  crt.bayes.control <- do.call("crt.bayes.control", crt.bayes.control)


  ## if only one point design was requested, then the weight can only be one
  if (k == 1)
    ICA.control$equal_weight <- TRUE
  if (length(lx) != 1 && ICA.control$sym)
    stop("currently symetric property only can be applied to models with one variable")

  # for the constraints
  if (!is.null(const$ui)){
    if (nrow(const$ui) != length(const$ci))
      stop("length of 'ci' must be equal to the number of rows in 'ui'")
    if (is.null(const$coef))
      stop("'coef' is missing")
    const$use <- TRUE
  }else
    const$use = FALSE
  npred <- length(lx)
  if (is.null(npar))
    npar <- prior$npar

  temp_crt <- create_criterion_bayes(FIM =  fimfunc2, type = type[1], prior= prior, const = const, compound = compound, multiple = multiple,
                                     localdes = NULL, method = crt_method[1], npar = npar, crt.bayes.control = crt.bayes.control, IRTpars = IRTpars)


  truncated_standard <- hcubature(f = function(param) prior$fn(t(param)),
                                  lowerLimit = prior$lower,
                                  upperLimit =  prior$upper, vectorInterface = TRUE)$integral

  temp_psi <- create_Psi_bayes(type = type[1], prior = prior, FIM = fimfunc2, lp = prior$lower, up = prior$upper, npar = npar, truncated_standard = truncated_standard, const = const, sens.bayes.control = sens.bayes.control, compound = compound,
                               method = sens_method[1])
  #############################################################################*
  # # construct_penalty
  # if (const$use){
  #   const$pen_func <- construct_pen(const = const, npred = npred, k = k)
  #   ## requird
  #   const$pen_func_1point <- construct_pen(const = const, npred = npred, k = 1)
  # }
  # #############################################################################*


  temp3 <- return_ld_ud (sym = ICA.control$sym, equal_weight = ICA.control$equal_weight, k = k, npred = npred, lx = lx, ux = ux)
  initial <- check_initial(initial = initial, ld = temp3$ld, ud = temp3$ud)

  #############################################################################*
  ## making the arg list
  ## the variables that will be added to control not by user, but by mica
  arg <- list(lx = lx, ux = ux, k = k, ld = temp3$ld, ud = temp3$ud, type = type[1],
              npar = npar,
              initial = initial, ICA.control = ICA.control,
              sens.bayes.control = sens.bayes.control,
              crt.bayes.control = crt.bayes.control,
              FIM = fimfunc2,
              crfunc = temp_crt$crfunc,
              prior = prior, const = const,
              Psi_funcs = temp_psi,
              plot_3d = plot_3d[1],
              compound = compound,
              truncated_standard = truncated_standard,
              crt_method = crt_method,
              sens_method = sens_method)






  ## because of Vector Interface
  #############################################################################*


  ICA_object <- list(arg = arg,
                     evol = NULL)
  class(ICA_object) <- c("list", "bayes")
  #cat("bayes_inner ", get(".Random.seed")[2], "\n")
  out <- iterate.bayes(object = ICA_object, iter = iter)
  out$arg$time <- proc.time() -time1
  return(out)

}
######################################################################################################*
######################################################################################################*

# ################# Quadrature Functions
# RSquadrature <- function(p, mu, Sigma, Nr, Nq) {
#
#   if(!identical(length(mu), as.integer(p))) stop("length(mu) must equal p")
#   if(!identical(length(Sigma), as.integer(p * p))) stop("Sigma must have size p * p")
#   #generate orthogonal matrices
#   #P<- p+1
#   P<-p
#   Q <- array(dim=c(Nq,P,P))
#
#   for (i in 1:Nq) {
#     ran.mat<-matrix(rnorm(P*P),ncol=P)
#     qr <- qr(ran.mat)
#     Q[i,,]<-qr.Q(qr)
#   }
#
#
#   # compute L, Cholesky root of Sigma
#   # L <- sqrt(Sigma)
#   L <- t(chol(Sigma))
#
#   # compute simplex + weights
#   simp <-simplex(P)
#   simplex <- simp$simplex
#   w.s <- simp$w.s
#
#   #print(simplex)
#   #print(w.s)
#   # compute radial abscissae, store in vector r
#
#   r <- gaulag(p=P,Nr=Nr,its=1e6,precision=1e-6)
#   w.R <- cassity.weight(r,P)/gamma((P)/2)
#   r <- 2*r
#   #print(r)
#
#   #create arrays to store abscissae and weights, and put values for the zero abscissa
#   a<- matrix(mu,ncol=P)
#   #a[1,P] <<- exp(mu[P])
#
#   w.a <- NULL
#   w.a <- w.R[1]
#
#   # calculate remaining abscissae and weights
#   for (i in 1:Nr) {
#     for (j in 1:((P+1)*(P+2))) {
#       for (k in 1:Nq) {
#         # compute the abscissa in beta-log(sigma^2) space
#         theta <- mu + (r[i+1]^0.5)*L%*%Q[k,,]%*%simplex[j,]
#
#         # exponentiate log(sigma^2) to put on beta-sigma^2 scale
#         #ONLY FOR GLMM case
#         #theta[P] <- exp(theta[P])
#         a <- rbind(a,t(theta))
#         w.a <- c(w.a,w.R[i+1]*w.s[j]/Nq)
#       }
#     }
#   }
#   return(list(a = a, w = w.a))
# }
#
# ###############################################################################################################*
# ###############################################################################################################*
# RSquadrature.uniform <- function(p, limits, Nr, Nq) {
#
#   #generate orthogonal matrices
#   #P<- p+1
#   P<-p
#   Q <- array(dim=c(Nq,P,P))
#
#   for (i in 1:Nq) {
#     ran.mat<-matrix(rnorm(P*P),ncol=P)
#     qr <- qr(ran.mat)
#     Q[i,,]<-qr.Q(qr)
#   }
#
#
#
#   # compute simplex + weights
#   simp <-simplex(P)
#   simplex <- simp$simplex
#   w.s <- simp$w.s
#
#   #print(simplex)
#   #print(w.s)
#
#   # compute radial abscissae, store in vector r
#
#   r <- gaulag(p=P,Nr=Nr,its=1e6,precision=1e-6)
#   w.R <- cassity.weight(r,P)/gamma((P)/2)
#   r <- 2*r
#   #print(r)
#
#   d1<-limits[,2]-limits[,1]
#   d2<-limits[,1]
#
#   #create arrays to store abscissae and weights, and put values for the zero abscissa
#   a <- matrix(d1*pnorm(0)+d2,ncol=P)
#   #a[1,P] <<- exp(mu[P])
#
#   w.a <- NULL
#   w.a <- w.R[1]
#
#   # calculate remaining abscissae and weights
#   for (i in 1:Nr) {
#     for (j in 1:((P+1)*(P+2))) {
#       for (k in 1:Nq) {
#         # compute the abscissa in beta-log(sigma^2) space
#         theta <- d1*pnorm((r[i+1]^0.5)*Q[k,,]%*%simplex[j,])+d2
#
#         # exponentiate log(sigma^2) to put on beta-sigma^2 scale
#         #ONLY FOR GLMM case
#         #theta[P] <- exp(theta[P])
#         a <- rbind(a,t(theta))
#         w.a <- c(w.a,w.R[i+1]*w.s[j]/Nq)
#       }
#     }
#   }
#   return(list(a = a, w = w.a))
# }
# ###############################################################################################################*
# ###############################################################################################################*
# simplex <- function(p) {
#   V <- matrix(ncol=p,nrow=p+1)
#   for (i in 1:(p+1))
#   {
#     for (j in 1:p) {
#       if (j<i) { V[i,j] = -((p+1)/(p*(p-j+2)*(p-j+1)))^0.5 }
#       if (j==i) { V[i,j] = ((p+1)*(p-i+1)/(p*(p-i+2)))^0.5}
#       if (j>i) { V[i,j] = 0 }
#     }
#   }
#
#   # Form midpoints and project onto the sphere
#
#   midpoints <- matrix(ncol=p,nrow=p*(p+1)/2)
#   k<-1
#   for (i in 1:p) {
#     for (j in (i+1):(p+1)) {
#       midpoints[k,] = 0.5*(V[i,]+V[j,])
#       k<-k+1
#     }
#   }
#   proj.pts <- midpoints # gets correct dimensions
#   for (k in 1:(p*(p+1)/2))
#   {
#     norm <- (sum(midpoints[k,]^2))^0.5
#     proj.pts[k,] <- midpoints[k,]/norm
#     if(identical(proj.pts[k, ], NaN)) proj.pts[k,] <- 0
#   }
#
#   # Form extended simplex by adding in negative images of these points
#   simplex <- rbind(V,-V,proj.pts,-proj.pts)
#
#   # Compute the simplex weights
#   w.s <- vector(length=(p+1)*(p+2))
#   w.s[1:(2*(p+1))] <- p*(7-p)/(2*(p+1)^2*(p+2))
#   w.s[-(1:(2*(p+1)))] <- 2*((p-1)^2)/(p*(p+1)^2*(p+2))
#   return(list(simplex=simplex,w.s=w.s))
# }
#
# ###############################################################################################################*
# ###############################################################################################################*
# # John's Laguerre poly root finder
# # translated from the C code in Press et al.
# # I have checked this against the latter, TW 20.1.2010.
# gaulag<-function(p, Nr, its,precision)
# {
#   alpha<-p/2
#   a<-vector()
#   w<-vector()
#   for(i in 1:Nr){
#     if(i==1){z<-(1+alpha)*(3+0.92*alpha)/(1+2.4*Nr+1.8*alpha)}
#     if(i==2){z<-z+(15+6.25*alpha)/(1+2.4*Nr+1.8*alpha)}
#     if(i>2){ai<-i-2}
#     if(i>2){z<-z+((1+2.55*ai)/(1.9*ai)+1.26*ai*alpha/(1+3.5*ai))*(z-a[i-2])/(1+0.3*alpha)}
#
#     for(its in 1:its){
#       p1<-1;p2<-0
#       for(j in 1:Nr){
#         p3<-p2;p2<-p1
#         p1<-((2*j-1+alpha-z)*p2-(j-1+alpha)*p3)/j}
#
#       pp<-(Nr*p1-(Nr+alpha)*p2)/z
#       z1<-z
#       z<-z1-p1/pp
#       if(abs(z-z1)< precision) break
#     }
#     a[i]<-z
#   }
#   a<-c(0,a)
#   return(a)
# }
# ###############################################################################################################*
# ###############################################################################################################*
# #Evaluates the laguerre polynomial with parameters n and s at the value a.
# laguerre<-function(a,n,s)
# {
#   laguerre.matrix<-matrix(nrow=length(a), ncol=n+1)
#   laguerre.vector<-vector(length=length(a))
#   #Loops up the recurrence relation
#   for(j in 1:length(a)){
#     for(i in 0:n){
#       laguerre.matrix[j,i+1]<-factorial(n)*choose(n+s, n-i)*(-a[j])^i/factorial(i)
#
#
#     }
#     laguerre.vector[j]<-sum(laguerre.matrix[j,])
#   }
#   return(laguerre.vector)
# }
# ###############################################################################################################*
# ###############################################################################################################*
# # Reproduces the weights formula given by Cassity(1965). But takes
# # the number of parameters as input to make it easier for the function user.
# cassity.weight<-function(a,p)
# {
#   n<-length(a);s<-p/2-1
#   constant<-gamma(n)*gamma(n+s)/(n+s)
#   weight<-constant/(laguerre(a,n-1,s)^2)
#   weight[1]<-weight[1]*(1+s)#as described by cassity et al. 'incorporate factor 1+s
#   # TW: Yeah, I find that sentence a bit confusing.
#   # but multiplying by 1+s here gives us the right answers
#   # (checking against the tables...)
#   return(weight)
# }
# ###############################################################################################################*
###############################################################################################################*

######################################################################################################*
######################################################################################################*
sensbayes_inner <- function(formula,
                            predvars, parvars,
                            family = gaussian(),
                            x, w,
                            lx, ux,
                            fimfunc = NULL,
                            prior = list(),
                            crt_method = c("cubature", "quadrature"),
                            sens_method = c("cubature", "quadrature"),
                            sens.bayes.control = list(),
                            crt.bayes.control = list(),
                            type = c("D", "DPA", "DPM", "multiple"),
                            plot_3d = c("lattice", "rgl"),
                            plot_sens = TRUE,
                            const = list(ui = NULL, ci = NULL, coef = NULL),
                            compound = list(prob = NULL, alpha = NULL),
                            varlist = list(),
                            calledfrom = c("sensfuncs", "iter", "plot"),
                            npar = NULL,
                            calculate_criterion = TRUE,
                            silent = FALSE,
                            calculate_sens = TRUE,
                            ...){
  time1 <- proc.time()

  if (calledfrom[1]  == "sensfuncs"){
    fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                         predvars = predvars, parvars = parvars,
                                         family = family, lx =lx, ux = ux,
                                         iter = 1, k = length(w),
                                         paramvectorized = TRUE, prior = prior)



    if(missing(formula)){
      # it should retunr a list of length one like the formula that returns a list
      fimfunc2 <- function(x, w, param)
        lapply(1:dim(param)[1], FUN = function(j)fimfunc(x = x, w = w, param = param[j, ]), ...)
      #fimfunc(x = x, w = w, param = param,...)
    } else{
      # if (length(prior$lower) != length(parvars))
      if (length(prior$lower) != fimfunc_formula$num_unknown_param)
        stop("length of 'prior$lower' is not equal to the number of unknown (not fixed) parameters")
      # fim_localdes <- fimfunc_formula$fimfunc_formula
      fimfunc2 <- fimfunc_formula$fimfunc_formula ## is vectorized with respect to the parameters
    }

    sens.bayes.control <- do.call("sens.bayes.control", sens.bayes.control)
    crt.bayes.control <- do.call("crt.bayes.control", crt.bayes.control)

    #############################################################################*
    #for the constraints
    if (!is.null(const$ui)){
      if (nrow(const$ui) != length(const$ci))
        stop("length of 'ci' must be equal to the number of rows in 'ui'")
      if (is.null(const$coef))
        stop("'coef' is missing")
      const$use <- TRUE
    }else
      const$use = FALSE
    #############################################################################*

    npred <- length(lx)
    if (is.null(npar))
      npar <- prior$npar

    temp_crt <- create_criterion_bayes(FIM =  fimfunc2, type = type[1], prior = prior, const = const, compound = compound, multiple = multiple,
                                       localdes = NULL, method = crt_method[1], npar = npar, crt.bayes.control = crt.bayes.control)




    ## the following si required for cheking the equivalence theorem
    truncated_standard <- hcubature(f = function(param) prior$fn(t(param)),
                                    lowerLimit = prior$lower,
                                    upperLimit =  prior$upper, vectorInterface = TRUE)$integral

    temp_psi <- create_Psi_bayes(type = type[1], prior = prior, FIM = fimfunc2, lp = prior$lower,
                                 up = prior$upper, npar = npar, truncated_standard = truncated_standard,
                                 const = const, sens.bayes.control = sens.bayes.control,
                                 compound = compound, method = sens_method[1])





    # if(length(lx) == 1)
    #   Psi_x_plot <-  temp_psi$Psi_x_bayes ## for PlotPsi_x

    # it is necessary to distniguish between Psi_x for plotiing and finding ELB becasue in plotting for models with two
    # explanatory variables the function should be defined as a function of x, y (x, y here are the ploints to be plotted)
    # if(length(lx) == 2)
    #   Psi_x_plot <- temp_psi$Psi_xy_bayes

    vertices_outer <- make_vertices(lower = lx, upper = ux)

    varlist <-list(npred =  length(lx),
                   crfunc = temp_crt$crfunc,
                   Psi_x_bayes  = temp_psi$Psi_x_bayes,
                   Psi_xy_bayes  = temp_psi$Psi_xy_bayes,
                   #Psi_x_minus = Psi_x_minus,
                   plot_3d = plot_3d[1],
                   vertices_outer = vertices_outer,
                   npar = npar)

  }
  if (calledfrom[1] ==  "iter")
    if (is.null(varlist) || is.null(npar) || calculate_criterion)
      stop("Bug: 'varlist' and 'npar' must be given when you call the sensitivity function from 'iter'\n'calculate_criterion' can not be TRUE because it is overdoing!")

  ##########################################################################*
  # find the maximum of derivative function
  # OptimalityCheck <- directL(fn = Psi_x_minus_bayes,
  #                            lower = lx, upper = ux,
  #                            x = x, w = w,
  #                            nl.info = FALSE,
  #                            control=list(xtol_rel=.Machine$double.eps,
  #                                         maxeval = 1000))
  ## change later
  #maxeval
  #const, Psi_x_bayes

  if (calculate_sens){
    if (is.null(sens.bayes.control$x0))
      x0 <- (lx + ux)/2 else
        x0 <- sens.bayes.control$x0

      if(!silent){
        #if (calledfrom == "sensfuncs")
        #cat("\n********************************************************************\n")
        cat("Please be patient! it may take very long for Bayesian designs....\nCalculating ELB..................................................\n")
      }
      OptimalityCheck <- nloptr::nloptr(x0= x0, eval_f = Psi_x_minus_bayes, lb = lx, ub = ux,
                                        opts = sens.bayes.control$optslist,
                                        x = x, w = w, const = const, Psi_x_bayes = varlist$Psi_x_bayes)

      ##sometimes the optimization can not detect maximum  in the bound, so here we add the cost values on the bound

      if (const$use){
        pen_val  <- apply(varlist$vertices_outer, 1, FUN = const$pen_func_1point)
        vertices_outer <- varlist$vertices_outer[which(pen_val == 0), , drop = FALSE]
      }

      if (nrow(varlist$vertices_outer) != 0){
        check_vertices <- find_on_points(fn = varlist$Psi_x_bayes,
                                         points = varlist$vertices_outer,
                                         x = x, w = w)

        vertices_val <- check_vertices$minima[, dim(check_vertices$minima)[2]]
        ## minus because for optimality check we minimize the minus sensitivity function
        max_deriv <- c(-OptimalityCheck$objective, vertices_val)
      }else
        max_deriv <- c(-OptimalityCheck$objective)
      max_deriv <- max(max_deriv)
      ##########################################################################*

      pointval <- find_on_points(fn = varlist$Psi_x_bayes,
                                 points = matrix(x, ncol = varlist$npred),
                                 x = x, w = w)

      # Efficiency lower bound
      ELB <- varlist$npar/(varlist$npar + max_deriv)
      if ((!silent))
        cat("    Maximum of the sensitivity function is ", max_deriv, "\n    Efficiency lower bound (ELB) is ", ELB, "\n")
      if (ELB < 0)
        warning("'ELB' can not be negative.\n1) Provide the number of model parameters by argument 'npar' when any of them are fixed. \n2) Increase the value of 'maxeval' in 'optslist'")
      if (length(lx) == 2)
        Psi_plot <- varlist$Psi_xy_bayes else
          Psi_plot <- varlist$Psi_x_bayes

      #Psi_x, const, plot_3d  = "lattice"
      if (plot_sens){
        if(!silent)
          cat("Plotting the sensitivity function................................\n")
        PlotPsi_x_bayes(lower = lx,
                        upper =  ux,
                        Psi_x = Psi_plot,
                        x = x, w = w, const = const,
                        plot_3d = plot_3d[1])
      }
      object <- list(type = type,
                     max_deriv = max_deriv,
                     ELB = ELB,
                     pointval = pointval)
  }
  if (calculate_criterion){
    if(!silent)
      cat("Evaluating the criterion.........................................\n")
    object$crtval <- varlist$crfunc(q = c(x, w), npred = varlist$npred)$val
  }
  class(object) <- c("list", "sensbayes")
  time2 <- proc.time() - time1
  if(!silent)
    cat("Verification is done in",time2[3], "seconds!", "\nAdjust the control parameters in 'sens.bayes.control' for higher speed\n")
  #, "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  object$time <- time2[3]
  return(object)
}


