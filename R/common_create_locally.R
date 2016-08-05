

create_locally <-  function(fimchar, lx, ux){
  # input: fimchar as a character.
  # output: compute the D-criterion given the locally optimal design

  # 'locally_fn' is created as a function of q, param, lx, ux and the Fisher information matrix
  #  locally_fn' is the objective function that is used in 'locally_det' to find the locally optimal design.
  # in 'locally_fn' argument 'q' is the point or weights of the design space that we should find numerically.
  # 'LB' and 'UB' in 'locally_det' should be adjusted according to 'q'.
  # 'locally_fn' return the -log(det(FIM)), but it still depends on q as an argument.
  # but 'locally_det' is only the function of 'param', so it finds the LDOD either by analytical solution or numerical optimization  on 'locally_fn'.
  # if analytical solution for LDOD is available, we do not need 'locally_fn'.


  if (fimchar == "FIM_exp_2par") {
    if (lx == 0 && ux == 1)
      fimchar <- "FIM_exp_2par_closed_form" else
        fimchar <- "FIM_exp_2par_numerical"
  }
  # Models that we only know the LDOD equally supported, but we dont have knowledge about the design points
  if (fimchar == "unknown" || fimchar == "FIM_logistic" ||
      fimchar == "FIM_power_logistic" || fimchar == "FIM_compartmental" ||
      fimchar == "FIM_exp_2par_numerical"){
    ## we know the weights, but we dont know the points. so 'q' is the points and 'LB' should be the lower bound for these points
    locally_fn <- function(q, param, auxiliary){
      x <- q # design points
      np <- length(x) # number of parameters
      w <- rep(1/np, np) # equally weighted designs
      det_fim <- det2(auxiliary$fimfunc(x=x, w = w, param = param), logarithm=FALSE)
      if (det_fim <= 0)
        cr <- 1e24 else
          cr <- -log(det_fim )
      return(cr)
    }
    locally_det <- function(param, auxiliary){
      LB  <-rep(auxiliary$lx, auxiliary$npar) # lowerbound
      UB <- rep(auxiliary$ux, auxiliary$npar) # upper bound
      locally_design <- gosolnp2(fun = locally_fn,
                                 LB = LB,
                                 UB = UB,
                                 n.sim = auxiliary$control$n.sim,
                                 n.restarts =  auxiliary$control$n.restarts,
                                 control = list(trace = auxiliary$control$trace),
                                 auxiliary = auxiliary,
                                 param = param)

      denominator <- exp(-locally_design$values[length(locally_design$values)])
      return(denominator)
    }
    use_gosolnp <- TRUE
  }

  # Models that the UPPER bound 'ux' is one of the support points of LDOD
  # equally weighted
  if (
    fimchar == "FIM_michaelis2" ||
    fimchar == "FIM_michaelis3" ||
    fimchar == "FIM_second_poi" ||
    fimchar == "FIM_poi2" ||
    fimchar == "FIM_sig_emax"
  ){

    ##"Sigmoid-EMAX" : is form the following paper
    ##Dette, Holger, Viatcheslav B. Melas, and Weng Kee Wong.
    #"Optimal design for goodness-of-fit of the Michaelis-Menten enzyme kinetic function.
    #" Journal of the American Statistical Association 100.472 (2005): 1370-1381.

    locally_fn <- function(q, param, auxiliary){
      x <- c(q, auxiliary$ux)
      np <- length(x)
      w <- rep(1/np, np)
      det_fim <- det2(auxiliary$fimfunc(x=x, w = w, param = param),
                      logarithm=FALSE)
      if (det_fim <= 0)
        cr <- 1e24 else
          cr <- -log(det_fim ) #+ 5000*(sum(w) - 1)^2

      return(cr)

    }
    locally_det <- function(param, auxiliary){
      LB  <-rep(auxiliary$lx, auxiliary$npar-1)
      UB <- rep(auxiliary$ux, auxiliary$npar-1)

      locally_design <- gosolnp2(fun = locally_fn,
                                 LB = LB,
                                 UB = UB,
                                 n.sim = auxiliary$control$n.sim,
                                 n.restarts = auxiliary$control$n.restarts,
                                 control = list(trace = auxiliary$control$trace),
                                 rseed = auxiliary$control$rseed,
                                 auxiliary = auxiliary,
                                 param = param)

      denominator <- exp(-locally_design$values[length(locally_design$values)])
      return(denominator)
    }
  }

  # Models that the LOWER bound 'lx' is one of the lsupport points of the LDOD
  ## equal weighted
  if (fimchar == "FIM_nbin"){
    locally_fn <- function(q, param, auxiliary){

      x <- c(auxiliary$lx,q)
      np <- length(x)
      w <- rep(1/np, np)
      det_fim <- det2(auxiliary$fimfunc(x=x, w = w, param = param),
                      logarithm=FALSE)
      if (det_fim <= 0)
        cr <- 1e24 else
          cr <- -log(det_fim )
      return(cr)
    }
    locally_det <- function(param, auxiliary){
      LB  <-rep(auxiliary$lx, auxiliary$npar-1)
      UB <- rep(auxiliary$ux, auxiliary$npar-1)

      locally_design <- gosolnp2(fun = locally_fn,
                                 LB = LB,
                                 UB = UB,
                                 n.sim = auxiliary$control$n.sim,
                                 n.restarts = auxiliary$control$n.restarts,
                                 control = list(trace = auxiliary$control$trace),
                                 rseed = auxiliary$control$rseed,
                                 auxiliary = auxiliary,
                                 param = param)

      denominator <- exp(-locally_design$values[length(locally_design$values)])
      return(denominator)
    }
    use_gosolnp <- TRUE
  }


  # Models that the LOWER bound 'lx' and Upper bound 'ux' are in the support points of LDOD
  if (fimchar == "FIM_weibull" || ##four parameters
      fimchar ==  "FIM_richards"  ## four parameters
  ){

    ##For richards and Weibull
    ##Dette, Holger, and Andrey Pepelyshev. "Efficient experimental designs for sigmoidal growth models." Journal of statistical planning and inference 138.1 (2008): 2-17.

    locally_fn <- function(q, param, auxiliary){

      x <- c(auxiliary$lx, q, auxiliary$ux)
      np <- length(x)
      w <- rep(1/np, np)
      det_fim <- det2(auxiliary$fimfunc(x=x, w = w, param = param),
                      logarithm=FALSE)
      if (det_fim <= 0)
        cr <- 1e24 else
          cr <- -log(det_fim ) #+ 5000*(sum(w) - 1)^2
      return(cr)
    }
    locally_det <- function(param, auxiliary){
      LB  <-rep(auxiliary$lx, auxiliary$npar-2)
      UB <- rep(auxiliary$ux, auxiliary$npar-2)
      locally_design <-  gosolnp2(fun = locally_fn,
                                  LB = LB,
                                  UB = UB,
                                  n.sim = auxiliary$control$n.sim,
                                  n.restarts =  auxiliary$control$n.restarts,
                                  control = list(trace = auxiliary$control$trace),
                                  rseed = auxiliary$control$rseed,
                                  auxiliary = auxiliary,
                                  param = param)
      denominator <- exp(-locally_design$values[length(locally_design$values)])
      return(denominator)
    }
    use_gosolnp <- TRUE
  }

  if ( fimchar == "FIM_sig_emax1" ||
       fimchar == "FIM_sig_emax3"
  ){
    ##Dette, Holger, and Andrey Pepelyshev. "Efficient experimental designs for sigmoidal growth models." Journal of statistical planning and inference 138.1 (2008): 2-17.
    locally_fn <- function(q, param, auxiliary){
      x <- c(q, auxiliary$ux)
      np <- length(x)
      w <- rep(1/np, np)
      det_fim <- det2(auxiliary$fimfunc(x=x, w = w, param = param),
                      logarithm=FALSE)

      if (det_fim <= 0)
        cr <- 1e24 else
          cr <- -log(det_fim ) #+ 5000*(sum(w) - 1)^2
      return(cr)
    }
    locally_det <- function(param, auxiliary){
      LB  <-rep(auxiliary$lx, auxiliary$npar-1)
      UB <- rep(auxiliary$ux, auxiliary$npar-1)
      locally_design <-  gosolnp2(fun = locally_fn,
                                  LB = LB,
                                  UB = UB,
                                  n.sim = auxiliary$control$n.sim,
                                  n.restarts =  auxiliary$control$n.restarts,
                                  control = list(trace = auxiliary$control$trace, tol = .Machine$double.xmin),
                                  rseed = auxiliary$control$rseed,
                                  auxiliary = auxiliary,
                                  param = param)
      denominator <- exp(-locally_design$values[length(locally_design$values)])
      return(denominator)
    }
  }

  if (fimchar == "FIM_poi1"){
    ## see Locally D- and c-optimal designs for Poisson
    ##and negative binomial regression models
    ##C. Rodriguez-Torreblanca 7 J. M. Rodriguez-Diaz
    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      denominator <- det2(auxiliary$fimfunc(x = c(auxiliary$lx, min(auxiliary$ux, auxiliary$lx + 2/param[2])), w = c(.5, .5), param = param), logarithm = FALSE)
      return(denominator)
    }
    use_gosolnp <- FALSE
  }

  if (fimchar == "FIM_emax_3par"){
    ## see Dette, H., Kiss, C., Bevanda, M. and Bretz, F. (2010), Optimal designs for the emax, log-linear and exponential models.
    ## Biometrika, 97 513-518.
    ## equation4
    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      xstar <- (auxiliary$ux*(auxiliary$lx+param[3])+auxiliary$lx*(auxiliary$ux+param[3]))/(auxiliary$lx + param[3] + auxiliary$ux + param[3])

      denominator <- det2(auxiliary$fimfunc(x = c(auxiliary$lx, xstar, auxiliary$ux) , w = rep(1/3, 3), param = param), logarithm = FALSE)
      return(denominator)
    }
    use_gosolnp <- FALSE
  }

  if (fimchar == "FIM_loglin"){
    ## see Dette, H., Kiss, C., Bevanda, M. and Bretz, F. (2010), Optimal designs for the emax, log-linear and exponential models.
    ## Biometrika, 97 513-518.
    ## equation5
    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      xstar <- (auxiliary$ux + param[3]) * (auxiliary$lx + param[3]) * (log(auxiliary$ux + param[3]) - log(auxiliary$lx + param[3]))/(auxiliary$ux - auxiliary$lx) - param[3]
      denominator <- det2(auxiliary$fimfunc(x = c(auxiliary$lx, xstar, auxiliary$ux) , w = rep(1/3, 3), param = param) , logarithm = FALSE)
      return(denominator)
    }
    use_gosolnp <- FALSE
  }
  if (fimchar == "FIM_exp_3par"){
    ## see Dette, H., Kiss, C., Bevanda, M. and Bretz, F. (2010), Optimal designs for the emax, log-linear and exponential models.
    ## Biometrika, 97 513-518.
    ## equation6
    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      xstar <- ((auxiliary$ux - param[3])*exp(auxiliary$ux/param[3]) - (auxiliary$lx - param[3])*exp(auxiliary$lx/param[3]))/
        (exp(auxiliary$ux/param[3]) - exp(auxiliary$lx/param[3]))
      denominator <- det2(auxiliary$fimfunc(x = c(auxiliary$lx, xstar, auxiliary$ux) , w = rep(1/3, 3), param = param), logarithm = FALSE)
      return(denominator)
    }
    use_gosolnp <- FALSE
  }

  if (fimchar == "FIM_michaelis"){
    ## see Dette, H., Kiss, C., Bevanda, M. and Bretz, F. (2010), Optimal designs for the emax, log-linear and exponential models.
    ## Biometrika, 97 513-518.
    ## equation6

    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      xstar <- param[2]  /(2*param[2]/auxiliary$ux + 1) ##equation 6
      denominator <- det2(auxiliary$fimfunc(x = c(xstar, auxiliary$ux) , w = c(.5, .5), param = param) ,logarithm = FALSE)
      ##equation 7
      #denominator <- 1/(64*param[2]^2*((1+param[2])^6))
      return(denominator)

    }
    use_gosolnp <- FALSE
  }
  if (fimchar == "FIM_exp_2par_closed_form"){
    #     //ON THE NUMBER OF SUPPORT POINTS OF MAXIMIN AND BAYESIAN OPTIMAL DESIGNS1
    #     // BRAESS AND  DETTE (2007)
    #     //a+exp(-bx)
    #     //example3.3
    #x=[0, 1] only!!!!
    locally_fn <- NA
    locally_det <- function(param, auxiliary){

      # warning: in the original paper the locally optimal design has only been found for lx = 0, but here by numerical results we
      # saw that it is true for all lx
      denominator <- det2(auxiliary$fimfunc(x = c(auxiliary$lx, 1/param[2]) , w = c(.5, .5), param = param) ,logarithm = FALSE)

      return(denominator)

    }
    use_gosolnp <- FALSE
  }
  if (fimchar == "FIM_comp_inhibition"){

    ###Optimum Design of Experiments for Enzyme Inhibition Kinetic Models
    #Journal of Biopharmaceutical Statistics
    #Barbara Bogacka a , Maciej Patan b , Patrick J. Johnson c , Kuresh Youdim c & Anthony C. Atkinson d
    #table one
    ##lx[1] is for S
    ##lx[2] for I
    #     param <- c(7.29753, 4.38594, 2.58218)
    #     lx = c(0, 0)
    #     ux = c(30 , 60)
    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      #first dimention is for S and the second one is for I.
      S_min <- auxiliary$lx[1]
      S_max <- auxiliary$ux[1]
      I_min <- auxiliary$lx[2]
      I_max <- auxiliary$ux[2]
      V <- param[1]
      K_m <- param[2]
      K_ic <- param[3]
      s2 <-  max(S_min, S_max*K_m*(K_ic + I_min))/(2*K_m*K_ic + 2*K_m*I_min+S_max*K_ic)
      s3 <- max(S_min, min(K_m*(K_ic + I_max)/K_ic, S_max))
      i3 <- min((2*K_m*I_min+S_max*K_ic+K_m*K_ic)/K_m, I_max)
      x <- c(S_max, s2, s3, I_min, I_min, i3)
      denominator <- det2(auxiliary$fimfunc(x = x , w = rep(1/3, 3), param = param) ,logarithm = FALSE)

      return(denominator)
    }
    use_gosolnp <- FALSE
  }

  if (fimchar == "FIM_uncomp_inhibition"){
    ###Optimum Design of Experiments for Enzyme Inhibition Kinetic Models
    #Journal of Biopharmaceutical Statistics
    #Barbara Bogacka a , Maciej Patan b , Patrick J. Johnson c , Kuresh Youdim c & Anthony C. Atkinson d
    #table one
    ##lx[1] is for S
    ##lx[2] for I

    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      #first dimention is for S and the second one is for I.
      S_min <- auxiliary$lx[1]
      S_max <- auxiliary$ux[1]
      I_min <- auxiliary$lx[2]
      I_max <- auxiliary$ux[2]
      V <- param[1]
      K_m <- param[2]
      K_iu <- param[3]
      s2 <- max(S_min, S_max*K_m*K_iu/(S_max*I_min+S_max*K_iu+2*K_m*K_iu))
      i3 <- min((2*S_max*I_min + S_max*K_iu+K_m*K_iu)/S_max, I_max)
      x <- c(S_max, s2, S_max, I_min, I_min, i3)
      denominator <- det2(auxiliary$fimfunc(x = x , w = rep(1/3, 3), param = param) ,logarithm = FALSE)

      return(denominator)

    }
    use_gosolnp <- FALSE
  }



  if (fimchar == "FIM_noncomp_inhibition"){

    ###Optimum Design of Experiments for Enzyme Inhibition Kinetic Models
    #Journal of Biopharmaceutical Statistics
    #Barbara Bogacka a , Maciej Patan b , Patrick J. Johnson c , Kuresh Youdim c & Anthony C. Atkinson d
    #table one
    ##lx[1] is for S
    ##lx[2] for I

    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      #first dimention is for S and the second one is for I.
      S_min <- auxiliary$lx[1]
      S_max <- auxiliary$ux[1]
      I_min <- auxiliary$lx[2]
      I_max <- auxiliary$ux[2]

      V <- param[1]
      K_m <- param[2]
      K_ic <- param[3]
      s2 <- max(S_min, S_max*K_m/(S_max+2*K_m))
      i3 <- min(K_ic+2*I_min, I_max)
      x <- c(S_max, s2, S_max, I_min, I_min, i3)

      denominator <- det2(auxiliary$fimfunc(x = x , w = rep(1/3, 3), param = param) ,logarithm = FALSE)

      return(denominator)

    }
    use_gosolnp <- FALSE
  }


  if (fimchar == "FIM_mixed_inhibition"){

    ###Optimum Design of Experiments for Enzyme Inhibition Kinetic Models
    #Journal of Biopharmaceutical Statistics
    #Barbara Bogacka a , Maciej Patan b , Patrick J. Johnson c , Kuresh Youdim c & Anthony C. Atkinson d
    #table one
    ##lx[1] is for S
    ##lx[2] for I

    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      #first dimention is for S and the second one is for I.
      S_min <- auxiliary$lx[1]
      S_max <- auxiliary$ux[1]
      I_min <- auxiliary$lx[2]
      I_max <- auxiliary$ux[2]
      V <- param[1]
      K_m <- param[2]
      K_ic <- param[3]
      K_iu <- param[4]
      s2 <- max(S_min, S_max*K_m*K_iu*(K_ic+I_min)/(S_max*K_ic*I_min+S_max*K_ic*K_iu+2*K_m*K_iu*I_min+2*K_m*K_iu*K_ic))
      i3 <- min((2*S_max*K_ic*I_min + S_max*K_ic*K_iu+2*K_m*K_iu*I_min+K_m*K_iu*K_ic)/(K_m*K_iu+S_max*K_ic), I_max)
      i4 <- min(I_min + (sqrt((K_ic+I_min)*(K_m*K_ic*K_iu+K_m*K_iu*I_min+S_max*K_ic*K_iu+S_max*K_ic*I_min)/
                                (K_m*K_iu+S_max*K_ic))),  I_max )

      s4 <- max(-K_m*K_iu*(K_ic+2*I_min-i4)/(K_ic*(K_iu+2*I_min-i4)),
                S_min)


      x <- c(S_max, s2, S_max, s4, I_min, I_min, i3, i4)
      denominator <- det2(auxiliary$fimfunc(x = x , w = rep(1/4, 4), param = param) ,logarithm = FALSE)

      return(denominator)

    }
    use_gosolnp <- FALSE
  }

  if (fimchar == "FIM_exp_2par_censor1" || fimchar == "FIM_exp_2par_censor2"){
    # theorem 2 optimal design for two-parameter nonlinear models with application to survival model
    locally_fn <- function(q, param, auxiliary){
      if (param[2] < 0)
        x <- c(auxiliary$lx, q)
      if (param[2]>0)
        x <- c(q, auxiliary$ux)

      np <- length(x)
      w <- rep(1/np, np)
      ##... is get from global environment
      det_fim <- det2(auxiliary$fimfunc(x=x, w = w, param = param),
                      logarithm=FALSE)
      if (det_fim <= 0)
        cr <- 1e24 else
          cr <- -log(det_fim ) #+ 5000*(sum(w) - 1)^2
      return(cr)
    }
    locally_det <- function(param, auxiliary){
      LB  <-rep(auxiliary$lx, auxiliary$npar-1)
      UB <- rep(auxiliary$ux, auxiliary$npar-1)
      locally_design <- gosolnp2(fun = locally_fn,
                                 LB = LB,
                                 UB = UB,
                                 n.sim = auxiliary$control$n.sim,
                                 n.restarts = auxiliary$control$n.restarts,
                                 control = list(trace = auxiliary$control$trace),
                                 rseed = auxiliary$control$rseed,
                                 auxiliary = auxiliary,
                                 param = param)

      denominator <- exp(-locally_design$values[length(locally_design$values)])
      return(denominator)
    }
    use_gosolnp <- TRUE

  }
  if (fimchar == "FIM_exp_3par_censor1" || fimchar == "FIM_exp_3par_censor2"){
    # theorem 5.1 on optimal design for censored data schmidt and schwabe 2015
    locally_fn <- function(q, param, auxiliary){

      beta0 <- param[1]
      beta1 <- param[2]
      beta2 <- param[3]
      if (beta2 > 0){
        if ((-beta1/beta2) > auxiliary$lx + auxiliary$ux)
          x <- c(auxiliary$lx, q)
        if ((-beta1/beta2) <= auxiliary$lx + auxiliary$ux)
          x <- c(q, auxiliary$ux)
        if ((-beta1/beta2) == auxiliary$lx + auxiliary$ux)
          x <- c(auxiliary$lx, q, auxiliary$ux)
      }

      if (beta2 < 0){
        if ((-beta1/beta2) <= auxiliary$lx )
          x <- c(auxiliary$lx, q) else
            if ((-beta1/beta2) >=  auxiliary$ux)
              x <- c(q, auxiliary$ux) else
                x <- q

      }
      np <- length(x)
      w <- rep(1/np, np)
      ##... is get from global environment
      det_fim <- det2(auxiliary$fimfunc(x=x, w = w, param = param),
                      logarithm=FALSE)
      if (det_fim <= 0)
        cr <- 1e24 else
          cr <- -log(det_fim ) #+ 5000*(sum(w) - 1)^2
      return(cr)
    }
    locally_det <- function(param, auxiliary){
      # fisrt we make the lower bound and upper bound for region of unceratinty
      beta0 <- param[1]
      beta1 <- param[2]
      beta2 <- param[3]

      if (beta2 > 0){
        if ((-beta1/beta2) == auxiliary$lx + auxiliary$ux){
          #because we only need to find one point
          LB  <-rep(auxiliary$lx, auxiliary$npar-2)
          UB <- rep(auxiliary$ux, auxiliary$npar-2)
        }else{
          LB  <-rep(auxiliary$lx, auxiliary$npar-1)
          UB <- rep(auxiliary$ux, auxiliary$npar-1)
        }
      }

      if (beta2 < 0){
        if ((-beta1/beta2) <= auxiliary$lx  || (-beta1/beta2) >=  auxiliary$ux){
          LB  <-rep(auxiliary$lx, auxiliary$npar-1)
          UB <- rep(auxiliary$ux, auxiliary$npar-1)

        }else{
          #browser()
          LB  <-rep(auxiliary$lx, auxiliary$npar)
          UB <- rep(auxiliary$ux, auxiliary$npar)
        }
      }
      locally_design <- gosolnp2(fun = locally_fn,
                                 LB = LB,
                                 UB = UB,
                                 n.sim = auxiliary$control$n.sim,
                                 n.restarts = auxiliary$control$n.restarts,
                                 control = list(trace = auxiliary$control$trace),
                                 rseed = auxiliary$control$rseed,
                                 auxiliary = auxiliary,
                                 param = param)

      denominator <- exp(-locally_design$values[length(locally_design$values)])
      return(denominator)
    }
    use_gosolnp <- FALSE
  }

  if (fimchar == "FIM_logistic_1par"){
    locally_fn <- NA
    locally_det<- function(param, auxiliary){
      ## See Grasshof, Holling and Schwabe 2012, optimal designs for the Rasch model
      denominator <- 1/4
      return(denominator)


    }
    use_gosolnp <- FALSE

  }

  #return(list(locally_det = locally_det, locally_fn = locally_fn))
  return(list(locally_det = locally_det, use_gosolnp = use_gosolnp))
}





