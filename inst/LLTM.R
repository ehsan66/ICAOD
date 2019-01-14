######################################################################################################*
######################################################################################################*
# next project
LLTM_locally <- function(Q,
                         lx, ux,  iter, k,
                         inipars,
                         ICA.control = list(),
                         sens.minimax.control = list(),
                         initial = NULL,
                         npar = length(inipars),
                         plot_3d = c("lattice", "rgl"),
                         normalization = NULL){

  if (is.null(normalization))
    normalization <- -mean(apply(Q, 1, function(row1) sum(row1 * inipars)))
  if (length(inipars) != dim(Q)[2])
    stop("length of 'inipars must be equal to the number of columns of 'Q' matrix")

  Q <- cbind(1, Q)
  inipars <- c(normalization, inipars)
  npar <- npar +1
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

    if (length(param)==1)
      lmat <- matrix(sum(simplify2array(lmat)), nrow = 1) else
        lmat <- apply(simplify2array(lmat), c(1, 2), sum)

    return(lmat)
  }

  out <- minimax_inner(fimfunc = fim_LLTM,
                       lx = lx, ux = ux, lp = inipars, up = inipars, iter = iter, k = k,
                       ICA.control = ICA.control,
                       sens.minimax.control = sens.minimax.control,
                       crt.minimax.control = list(inner_space = "locally"),
                       type = "locally",
                       initial = initial,
                       localdes = NULL,
                       npar = npar,
                       robpars = list(),
                       crt_type = "D",
                       multipars = list(),
                       plot_3d = plot_3d[1])

  return(out)
}







######################################################################################################*
######################################################################################################*
# next project
LLTM_minimax <- function(Q,
                         lx, ux,  iter, k,
                         lp, up,
                         n.grid = 0,
                         ICA.control = list(),
                         crt.minimax.control = list(),
                         sens.minimax.control = list(),
                         initial = NULL,
                         npar = length(inipars),
                         plot_3d = c("lattice", "rgl"),
                         normalization = NULL){



  if (!is.numeric(n.grid) || n.grid < 0)
    stop("value of 'n.grid' must be >= 0")
  crt.minimax.control <- do.call("crt.minimax.control", crt.minimax.control)
  # check if the length of x0 be equal to the length of lp and up!!
  #minimax.control$inner_maxeval <-  param_maxeval

  if (n.grid>0){
    crt.minimax.control$param_set <- make_grid(lp = lp, up = up, n.grid = n.grid)
    crt.minimax.control$inner_space <- "discrete"
  }else
    crt.minimax.control$inner_space <- "continuous"

  Q <- cbind(1, Q)

  npar <- npar +1
  lp <- c(normalization, lp)
  up <- c(normalization, up)
  wlambda <- function(x, w, param, q){
    # x is the design points that are the values for the ability parameters
    #qparam <- apply(q * param, 2, sum)
    qparam <- -sum(q * param)
    # argument of lambda for each param vector
    arg_lambda <- sapply(1:length(qparam), function(k)sum(w * exp(qparam[k] + x)/(1 + exp(qparam[k] + x))^2))
  }

  fim_LLTM <- function(x, w, param){
    Qmat <- Q
    # param <- c(normalization, param)
    nitems <- nrow(Qmat)
    lmat <- lapply(1:nitems, FUN = function(j)wlambda(q = Qmat[j, ],x = x, w = w, param = param) * Qmat[j, ] %*% t(Qmat[j, ]))

    if (length(param)==1)
      lmat <- matrix(sum(simplify2array(lmat)), nrow = 1) else
        lmat <- apply(simplify2array(lmat), c(1, 2), sum)

    return(lmat)
  }
  out <- minimax_inner(fimfunc = fim_LLTM,
                       lx = lx, ux = ux, lp = lp, up = up, iter = iter, k = k,
                       ICA.control = ICA.control,
                       sens.minimax.control = sens.minimax.control,
                       crt.minimax.control = crt.minimax.control,
                       type = "minimax",
                       initial = initial,
                       localdes = NULL,
                       npar = npar,
                       robpars = list(),
                       crt_type = "D",
                       multipars = list(),
                       plot_3d = plot_3d[1])

  return(out)

  ######################################################################################################*
  ######################################################################################################*
  # @title  Bayesian D-Optimal Designs for Linear Logistic Test Model (LLTM)
  # @description
  #  Finds (pseudo) Bayesian D-optimal designs for the LLTM.
  #  It should be used when the user assumes a (truncated) prior distribution for the unknown basic parameters.
  #
  # @inheritParams bayes
  # @inherit bayes return
  # @param Q Matrix of weights for basic rules.
  # @export
  # next project
  LLTM_bayes <- function(prior,
                         Q = NULL,
                         lx,
                         ux,
                         iter,
                         k,
                         ICA.control =  list(),
                         crt.bayes.control = list(),
                         sens.bayes.control = list(),
                         crt_method = c("cubature", "quadrature"),
                         sens_method = c("cubature",  "quadrature"),
                         initial = NULL,
                         npar = NULL,
                         plot_3d = c("lattice", "rgl")
                         ,normalization = NULL) {

    npar <- prior$npar + 1
    if (!is.numeric(npar))
      stop("'npar' is the number of parameters and must be numeric")

    Q <- cbind(1, Q)
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
      if (is.null(normalization))
        normalization <- -mean(apply(Qmat[, -1], 1, function(row1) sum(row1 * inipars)))
      param <- c(normalization, param)
      lmat <- lapply(1:nitems, FUN = function(j)wlambda(q = Qmat[j, ],x = x, w = w, param = param) * Qmat[j, ] %*% t(Qmat[j, ]))
      lmat <- apply(simplify2array(lmat), c(1, 2), sum)
      return(lmat)
    }
    output <-  bayes_inner(fimfunc = fim_LLTM,
                           lx = lx,
                           ux = ux,
                           type = "D",
                           iter = iter,
                           k = k,
                           npar = npar,
                           prior = prior,
                           compound = list(prob = NULL, alpha = NULL),
                           multiple.control = list(),
                           crt_method = crt_method[1],
                           sens_method = sens_method[1],
                           ICA.control =  ICA.control,
                           crt.bayes.control = crt.bayes.control,
                           sens.bayes.control = sens.bayes.control,
                           initial = initial,
                           plot_3d = plot_3d[1],
                           const = list(ui = NULL, ci = NULL, coef = NULL))


    return(output)
  }
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
    nw <- createNIGrid(dim = length(prior$lower), type = sens.bayes.control$quadrature$type,
                       level = sens.bayes.control$quadrature$level,
                       ndConstruction = sens.bayes.control$quadrature$ndConstruction)

    if (sens.bayes.control$quadrature$type == "GLe")
      rescale(nw, domain = matrix(c(prior$lower, prior$upper), ncol=2))
    if (sens.bayes.control$quadrature$type == "GHe")
      rescale(nw, m = prior$mu, C =  prior$sigma)

    if (type == "D"){
      Psi_x_integrand_D  <- function(x1, param,  FIM, x, w, prior_func){
        ## param: as matrix: EACH row is one set of parameters
        ## return a matrix with 1 * length(param) dimension
        FIM_x <- FIM(x = x, w = w, par = param)
        FIM_x1 <- FIM(x = x1, w = 1, par = param)
        deriv_integrand <- sapply(1:nrow(param), FUN = function(j)sum(diag(solve(FIM_x[[j]])%*%FIM_x1[[j]]))) * prior_func(param)
        if (dim(param)[1] != dim(deriv_integrand)[1])
          deriv_integrand  <- t(deriv_integrand)
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
        out <- quadrature(f = Psi_x_integrand_D, grid = nw, x=x, w=w, FIM = FIM, x1 = x1, prior_func = prior$fn)
        #return(out/truncated_standard  - npar)
        return(out  - npar)
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
        out <- quadrature(f = Psi_x_integrand_D, grid = nw, x=x, w=w, FIM = FIM,  x1 = c(x1, y1), prior_func = prior$fn)
        return(-(out-npar))
      }
    }
    if (type == "DPA" || type == "DPAM"){
      Psi_x_integrand_DP <- function(x1, param,  FIM,  x, w, prior_func){
        # NOTE: in cubature when integrand is vectorized:
        #  1- the input of the integrand (param) is a MATRIX with length(mean) or length(lp) number of rows and l number of columns  (l is for vectorization)
        #  2- the integrand should return a matrix of 1 * l (vectorization).
        ## param: as matrix: EACH COLUMN is one set of parameters

        if (type == "DPA"){
          FIM_x <- FIM(x = x, w = w, par = param)
          FIM_x1 <- FIM(x = x1, w = 1, par = param)
          if (compound$alpha != 0){
            deriv_integrand <-  sapply(1:nrow(param), FUN = function(j)compound$alpha/npar * sum(diag(solve(FIM_x[[j]]) %*% FIM_x1[[j]])) +
                                         (1-compound$alpha) * (compound$prob(x1, param[j, ])- sum(w * compound$prob(x, param[j, ])))/sum(w * compound$prob(x, param[j, ]))) * prior_func(param)
          }else{
            # when alpha = 0
            deriv_integrand <-  sapply(1:nrow(param), FUN = function(j) (compound$prob(x1, param[j, ])- sum(w * compound$prob(x, param[j, ])))/sum(w * compound$prob(x, param[j, ]))) #* prior_func(param)
          }
        }else
          if (type == "DPM"){
            stop("Bug:not completed for DPM")
            deriv_integrand <- apply(param, 2, FUN = function(col_par)compound$alpha/npar * sum(diag(solve(FIM(x = x, w = w, par = col_par)) %*% FIM(x = x1, w = 1, par = col_par))) + (1-compound$alpha) * (compound$prob(x1, col_par)- min(compound$prob(x, col_par)))/min(compound$prob(x, col_par))) * prior_func(t(param))
          }else
            stop("BUG: check 'type'")
        if (dim(param)[1] != dim(deriv_integrand)[1])
          deriv_integrand  <- t(deriv_integrand)
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
        out <- quadrature(f = Psi_x_integrand_DP, grid = nw, x=x, w=w, FIM = FIM, x1 = x1, prior_func = prior$fn)
        return(out - compound$alpha)
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
        out <- quadrature(f = Psi_x_integrand_DP, grid = nw, x=x, w=w, FIM = FIM, x1 = c(x1, y1), prior_func = prior$fn)
        return(-(out - compound$alpha))
      }
    }


  }
  return(list(Psi_xy_bayes = Psi_xy_bayes, Psi_x_bayes = Psi_x_bayes))
}
######################################################################################################*
######################################################################################################*
#############################################################################################################*
#############################################################################################################*
create_criterion_bayes <- function(FIM, type, prior, compound, const, multiple, localdes = NULL, method, npar, crt.bayes.control, IRTpars, is.only.w, only_w_varlist = NULL){
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


    if(type == "multiple")
      stop("BUG: No Bayesian multiple objective optimal designs is implemented yet!")

    cr_integrand <- switch(type, "D" = cr_integrand_D, "DPA" = cr_integrand_DPA, "DPM" = cr_integrand_DPM, "D_LLTM" = crt_LLTM)
    #cr_integrand <- switch(type, "D" = cr_integrand_D, "DPA" = cr_integrand_DPA, "DPM" = cr_integrand_DPM, "multiple" = cr_integrand_multiple)
    crfunc_bayesian  <- function(q, npred) {
      if (!is.only.w){
        lq <- length(q) # q is the design points and weights
        pieces <- lq / (npred + 1)
        x_ind <- 1:(npred * pieces)
        w_ind <- (x_ind[length(x_ind)] + 1):lq
        x <- q[x_ind]
        w <- q[w_ind]
      }else{
        w <- q
        x <- only_w_varlist$x
      }
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
    if (type == "D")
      cr_integrand_D<- function(param, x, w){
        bcrfunc1 <- -1 * sapply(FIM(x = x, w = w, param = param), det2, logarithm= TRUE) * prior_func(param)
        if (dim(param)[1] != dim(bcrfunc1)[1])
          bcrfunc1 <- t(bcrfunc1)
        return(bcrfunc1)
      }

    if (type == "DPA")
      cr_integrand_DPA <- function(param, x, w){
        bcrfunc1 <- -(compound$alpha/4  * sapply(FIM(x = x, w = w, param = param), det2, logarithm= TRUE) +
                        apply(param, 1, function(row_par)(1 -compound$alpha) * log(sum(w * compound$prob(x = x, param = row_par))))) * prior_func(param)
        if(any(is.infinite(bcrfunc1)))
          bcrfunc1[(is.infinite(bcrfunc1))] <- 0
        if (dim(param)[1] != dim(bcrfunc1)[1])
          bcrfunc1 <- t(bcrfunc1)
        return(bcrfunc1)
      }

    cr_integrand <- switch(type, "D" = cr_integrand_D, "DPA" = cr_integrand_DPA, "DPM" = cr_integrand_DPM, "D_LLTM" = crt_LLTM)

    nw <- createNIGrid(dim = length(prior$lower), type = crt.bayes.control$quadrature$type,
                       level = crt.bayes.control$quadrature$level,
                       ndConstruction = crt.bayes.control$quadrature$ndConstruction)
    if (crt.bayes.control$quadrature$type == "GLe")
      rescale(nw, domain = matrix(c(prior$lower, prior$upper), ncol=2))
    if (crt.bayes.control$quadrature$type == "GHe")
      rescale(nw, m = prior$mu, C =  prior$sigma)

    crfunc_bayesian  <- function(q, npred) {
      if (!is.only.w){
        lq <- length(q) # q is the design points and weights
        pieces <- lq / (npred + 1)
        x_ind <- 1:(npred * pieces)
        w_ind <- (x_ind[length(x_ind)] + 1):lq
        x <- q[x_ind]
        w <- q[w_ind]
      }else{
        w <- q
        x <- only_w_varlist$x
      }
      out <- quadrature(f = cr_integrand, grid = nw, x=x, w=w)
      if (out == Inf || is.nan(out) || out == -Inf )
        out <- 1e+20
      if (const$use)
        val <- out + 5000 * (sum(w) - 1)^2 + const$pen_func(x) else
          val <- out + 5000 * (sum(w) - 1)^2
      return(list(val = val, fneval = dim(nw$nodes)[1]))
    }
    crfunc <- crfunc_bayesian
  }

  #####################################################################*
  return(list(crfunc = crfunc))
}
###############################################################################################################*
###############################################################################################################*
###############################################################################################################*
###############################################################################################################*
plotLLTM_1point_bayes <- function(prior, Q, lx, ux, optimalx){


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



######################################################################################################*
######################################################################################################*
#'@title Calculates Relative Efficiency for Locally Optimal Designs
#' @description
#' Given a vector of initial estimates for the parameters, this function calculates the D-and PA- efficiency of a design \eqn{\xi_1} with respect to a design \eqn{\xi_2}.
#' Usually, \eqn{\xi_2} is an  optimal design.
#'
#'
#' @details For a known \eqn{\theta_0}, relative D-efficiency is
#' \deqn{exp(\frac{log|M(\xi, \theta_0)| - log|M(\xi^*, \theta_0)|}{npar})}{
#' exp\{(log|M(\xi, \theta_0)| - log|M(\xi*, \theta_0)|)/npar\}.}
#' The relative P-efficiency is
#' \deqn{\exp(\log(\sum_{i=1}^k w_ip(x_i, \theta_0) - \log(\sum_{i=1}^k w^*_ip(x^*_i, \theta_0))}{
#' exp(log (\sum w_i p(x_i, \theta_0) - log(\sum w*_i p(x*_i, \theta_0)),
#' }
#' where \eqn{x^*}{x*} and \eqn{w^*}{w*} are the support points and the corresponding weights of the optimal design, respectively.
#'
#'  The argument  \code{x} is the vector of design points.
#'  For design points with more than one dimension,
#'    it is a concatenation of the design points, but \strong{dimension-wise}.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  (three-dimensional) design points of a two-point optimal design are
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     Then, the argument \code{x} is equivalent to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'
#' @export
#' @inheritParams locallycomp
#' @param x Vector of design (support) points of \eqn{\xi_1}. See 'Details' of \code{\link{leff}}.
#' @param w Vector of corresponding design weights for \code{x}.
#' @param xopt Vector of design (support) points of optimal design (\eqn{\xi_2}). Similar to \code{x}.
#' @param wopt Vector of corresponding design weights for \code{xopt}.
#' @param type A character. \code{"D"} denotes D-efficiency and \code{"PA"} denotes average P-efficiency.
#' @return A value between 0 and 1.
#' @references McGree, J. M., Eccleston, J. A., and Duffull, S. B. (2008). Compound optimal design criteria for nonlinear models. Journal of Biopharmaceutical Statistics, 18(4), 646-661.
#' @example inst/examples/leff_examples.R
leff <- function(formula,
                 predvars,
                 parvars,
                 family = gaussian(),
                 inipars,
                 type = c("D", "PA"),
                 fimfunc = NULL,
                 xopt, wopt, x, w,
                 npar = length(inipars),
                 prob = NULL,
                 Q = NULL,
                 normalization = NULL){
  ## bayesian D-efficiency
  ### relative efficieny of x with respect to xopt
  ## if the user use small p

  if (!(type[1] %in% c("D", "PA")))
    stop("'type' must be 'D' or 'PA'")
  npred <- length(x)/length(w)
  if (!is.null(Q)){
    if (is.null(normalization))
      normalization <- -mean(apply(Q, 1, function(row1) sum(row1 * inipars)))
    Q <- cbind(1, Q)
    fimfunc <- create_FIM_LLTM(Q = Q, normalization = normalization)
    # wlambda <- function(x, w, param, q){
    #   # x is the design points that are the values for the ability parameters
    #   #qparam <- apply(q * param, 2, sum)
    #   qparam <- -sum(q * param)
    #   # argument of lambda for each param vector
    #   arg_lambda <- sapply(1:length(qparam), function(k)sum(w * exp(qparam[k] + x)/(1 + exp(qparam[k] + x))^2))
    # }
    #
    # fim_LLTM <- function(x, w, param){
    #   Qmat <- Q
    #   nitems <- nrow(Qmat)
    #   if (!is.null(normalization))
    #     param <- c(normalization, param)
    #   lmat <- lapply(1:nitems, FUN = function(j)wlambda(q = Qmat[j, ],x = x, w = w, param = param) * Qmat[j, ] %*% t(Qmat[j, ]))
    #   lmat <- apply(simplify2array(lmat), c(1, 2), sum)
    #   return(lmat)
    # }
    # fimfunc <- fim_LLTM
  }

  fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =rep(0, npred), ux = rep(1, npred), iter = 1, k = 1,
                                       paramvectorized = FALSE,
                                       prior = NULL, x = NULL)


  if (!missing(formula)){
    if (length(inipars) != length(parvars))
      stop("length of 'inipars' is not equal to the length of 'parvars'")
  }

  if(missing(formula)){
    # to handle ...
    fimfunc2 <- function(x, w, param)
      fimfunc(x = x, w = w, param = param)
    #fimfunc(x = x, w = w, param = param,...)
  } else{
    fimfunc2 <- fimfunc_formula$fimfunc_formula ## can be vectorized with respect to parameters!
  }
  if (type[1] == "D")
    releff <- (det2(fimfunc2(x = x, w = w, param = inipars))/
                 det2(fimfunc2(x = xopt, w = wopt, param = inipars)))^(1/npar)

  if (type[1] == "PA"){

    if (is.null(prob))
      stop("'prob' must be given for P-optimality")

    if (is.formula(prob)){
      prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
    }else{
      if (!is.function(prob))
        stop("'prob' must be either a function or a formula")
      if (!formalArgs(prob) %in% c("x", "param"))
        stop("arguments of 'prob' must be 'x' and 'param'")
    }
    releff <- sum(w * prob(x = x, param = inipars))/sum(wopt * prob(x = xopt, param = inipars))
  }


  return(releff)
}

######################################################################################################*
######################################################################################################*
#' @title  Calculates Relative Efficiency for Bayesian Optimal Designs
#' @description
#' Given a prior distribution for the parameters, this function calculates the Bayesian D-and PA- efficiency of a design \eqn{\xi_1} with respect to a design \eqn{\xi_2}.
#' Usually, \eqn{\xi_2} is an optimal design.
#' This function is especially useful for investigating the robustness of Bayesian optimal designs under different prior distributions (See 'Examples').
#'
#' @inheritParams leff
#' @inheritParams bayes
#' @details
#' See Masoudi et al. (2018) for formula details (the paper is under review and will be updated as soon as accepted).
#'
#'  The argument  \code{x} is the vector of design points.
#'  For design points with more than one dimension,
#'    it is a concatenation of the design points, but \strong{dimension-wise}.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  (three-dimensional) design points of a two-point optimal design are
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     Then, the argument \code{x} is equivalent to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'
#'
#' @export
#' @example inst/examples/beff_examples.R
beff <- function(formula,
                 predvars,
                 parvars,
                 family = gaussian(),
                 prior,
                 fimfunc = NULL,
                 xopt, wopt, x, w,
                 crt.bayes.control = list(),
                 npar = NULL,
                 type = c("D", "PA"),
                 prob = NULL,
                 Q = NULL,
                 normalization = NULL){
  ## bayesian D-efficiency
  ### relative efficieny of x with respect to xopt
  if (!(type[1] %in% c("D", "PA")))
    stop("'type' must be 'D' or 'PA'")
  if (is.null(npar)){
    if (!missing(formula))
      npar <- length(parvars)
    else
      npar <- prior$npar
  }
  if (!is.null(Q)){
    if (is.null(normalization))
      warning("you have requested LLTM. No normalization constant was given") else{
        Q <- cbind(1, Q)
      } # else
    fimfunc <- create_FIM_LLTM(Q = Q, normalization = normalization)
  }


  ## only to pass the check_common_eargs
  npred <- length(x)/length(w)
  fimfunc_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =rep(0, npred), ux = rep(1,npred),
                                       iter = 1, k = length(w),
                                       paramvectorized = FALSE, prior = prior,
                                       x = NULL)
  if(missing(formula)){
    # to handle ...
    fimfunc2 <- function(x, w, param){
      fimfunc(x = x, w = w, param = param)
      #fimfunc(x = x, w = w, param = param,...)
    }
  } else{

    #if (length(prior$lower) != length(parvars))
    if (length(prior$lower) != fimfunc_formula$num_unknown_param)
      stop("length of 'prior$lower' is not equal to the number of unknown (not fixed) parameters")
    # fim_localdes <- fimfunc_formula$fimfunc_formula
    fimfunc2 <- fimfunc_formula$fimfunc_formula ## can be vectorized with respect to parameters!
  }
  crt.bayes.control <- do.call("crt.bayes.control", crt.bayes.control)

  if (type[1] == "D")
    cr_integrand <- function(param, x, w){
      # bcrfunc1 <- apply(param, 2,
      #                   FUN = function(col_par)-det2(fimfunc2(x = x, w = w, param = col_par), logarithm = TRUE)) * prior$fn(t(param))
      bcrfunc1 <- apply(param, 2,
                        FUN = function(col_par)-det2(fimfunc2(x = x, w = w, param = col_par), logarithm = TRUE)) * prior$fn(t(param))

      dim(bcrfunc1) <- c(1, length(bcrfunc1))
      return(bcrfunc1)
    }
  if (type[1] == "PA"){
    if (is.formula(prob)){
      prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
    }else{
      if (!is.function(prob))
        stop("'prob' must be either a function or a formula when type is 'PA'")
      if (!formalArgs(prob) %in% c("x", "param"))
        stop("arguments of 'prob' must be 'x' and 'param'")
    }
    cr_integrand <- function(param, x, w){
      bcrfunc1 <- (apply(param, 2, function(col_par) -log(sum(w * prob(x = x, param = col_par))))) * prior$fn(t(param))
      dim(bcrfunc1) <- c(1, length(bcrfunc1))
      return(bcrfunc1)
    }

  }
  crfunc_bayesian  <- function(x, w, maxEval, tol) {

    out <- hcubature(f = cr_integrand, lowerLimit = prior$lower,
                     upperLimit = prior$upper,
                     vectorInterface = TRUE,
                     x = x, w = w, tol = tol, maxEval = maxEval)

    val <- out$integral
    return(list(val = val, fneval = out$functionEvaluations))
  }

  if (type[1] == "D")
    releff<- exp((crfunc_bayesian(x = xopt, w = wopt, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val-
                    crfunc_bayesian(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)/npar)
  # if (type[1] == "PA")
  #   releff<-  (crfunc_bayesian(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val/
  #                   crfunc_bayesian(x = xopt, w = wopt, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)
  if (type[1] == "PA")
    releff <- exp(crfunc_bayesian(x = xopt, w = wopt, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val - crfunc_bayesian(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val)

  # releff<- crfunc_bayesian_D(x = xopt, w = wopt, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val/
  #                crfunc_bayesian_D(x = x, w = w, maxEval = crt.bayes.control$cubature$maxEval, tol = crt.bayes.control$cubature$tol)$val

  return(releff)
}
######################################################################################################*
######################################################################################################*

######################################################################################################*
######################################################################################################*

######################################################################################################*
######################################################################################################*

