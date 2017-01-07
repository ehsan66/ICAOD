# roxygen
#' Update an object of class 'ICA'
#'
#' Function \code{iterate.ICA} takes an object of class \code{"ICA"} and runs the ICA algorithm and updates the results.
#'
#' @param object an object of class 'ICA'.
#' @param iter number of iterations.
#' @seealso \code{\link{mica}}
#' @importFrom nloptr directL
#' @export
iterate.ICA <- function(object, iter){

  if (all(class(object) != c("list", "ICA")))
    stop("''object' must be of class 'ICA'")
  if (missing(iter))
    stop("'iter' is missing")

  arg <- object$arg
  control <- object$arg$control
  evol <- object$evol
  n_independent <- length(arg$lx)
  lp <- arg$lp
  up <- arg$up

  ## all fo the types for optim_on_average will be set to be equal to "optim_on_average"
  ## but the arg$type remains unchanged to be used in equivalence function!!
  if( grepl("on_average", arg$type))
    type = "optim_on_average" else
      type <- arg$type
  # warning: no arg$type must be used further

  if (type == "optim_on_average")
    npar <- dim(arg$param)[2] else
      npar <- length(arg$lp)
  if (control$equal_weight)
    w_equal <- rep(1/arg$k, arg$k)

  ###############################################################
  ## multi_locally is the same as locally in update!

  if (type == "multiple_locally"){
    type <- "locally"
    ## rewuired for setting the title of plots
    multi_type <- TRUE
  }else
    multi_type <- FALSE

  # if (type == "multiple_minimax")
  #   type <- "minimax"

  if (!(type %in% c("minimax", "standardized", "locally", "optim_on_average")))
    stop("bug: the type must be 'minimax' or 'standardized' or 'locally' or 'optim_on_average' in 'iterate.ICA\nset  'multiple_locally' to 'locally'")
  # because they have the same configuration. But we need to know the multi becasue of the verifying and plot methods!
  ###################################################################################

  if (type == "locally")
    param_locally <- up
  if (type == "optim_on_average")
    param_set <- arg$param


  ##############################################################################
  ###finding if there is any fixed parameters.
  # only if type != "locally"
  #if (type != "locally" & control$inner_space != "vertices" & control$inner_space != "discrete"){
  if (type != "locally" && type != "optim_on_average"){
    # here we search if one of the parameters are fixed. then we pass it to the optimization function in the inner problem because otherwise it may casue an error.
    any_fixed <- sapply(1:length(lp), function(i) lp [i] == up[i]) # is a vector
    if (any(any_fixed)){
      is_fixed <- TRUE
      fixedpar_id <- which(any_fixed)
      fixedpar <- lp[fixedpar_id]
      lp <- lp[-fixedpar_id]
      up <- up[-fixedpar_id]
      ## warning: lp and up are channged here if 'any_fix == TRUE'
    }else{
      fixedpar <- NA
      fixedpar_id <- NA
      is_fixed <- FALSE
    }
  }else{
    fixedpar <- NA
    fixedpar_id <- NA
    is_fixed <- FALSE
  }
  if(control$inner_space == "discrete"){
    if(!is.na(fixedpar))
      discrete_set <- control$param_set[, -fixedpar_id, drop = FALSE] else
        discrete_set <- control$param_set
  }
  #############################################################################

  #############################################################################
  # plot setting
  #plot_cost <- control$plot_cost
  #plot_deriv <- control$plot_deriv
  legend_place <- "topright"
  legend_text <- c( "Best Imperialist", "Mean of Imperialists")
  line_col <- c("firebrick3", "blue4")
  if (type == "minimax")
    title1 <- "cost value"
  if (type == "standardized")
    title1 <- "minimum efficiency"
  if (type == "locally" || type == "optim_on_average")
    title1 <-  "log determinant of inverse of FIM"
  if (multi_type)
    title1 <- "criterion value"
  ################################################################################

  ## In last iteration the check functions should be applied??
  check_last <- ifelse(control$equivalence_every != FALSE, TRUE, FALSE)

  ################################################################################
  ### re-defimimg crfunc to handle fixed parameters.
  crfunc <- arg$crfunc
  if (is_fixed){
    crfunc2 <- function(param, q, fixedpar = NA, fixedpar_id = NA, n_independent){
      # if (any(!is.na(fixedpar))){
      #   if (any(is.na(fixedpar_id)))
      #     stop("'fixedpar' index is missing.")
      param_new <- rep(NA, length(param) + length(fixedpar))
      param_new[fixedpar_id] <- fixedpar
      param_new[-fixedpar_id] <- param
      param <- param_new
      #}
      out <- crfunc(param = param, q = q, n_independent = n_independent)
      return(out)
    }
  }else{
    crfunc2 <- function(param, q, fixedpar = NA, fixedpar_id = NA, n_independent){
      # no use for fixedpar  and fixedpar_id = NA
      out <- crfunc(param = param, q = q, n_independent = n_independent)
      return(out)
    }
  }
  ################################################################################

  ################################################################################
  ### Psi as a function of x and x, y for plotting. Psi_x defined as minus psi to find the minimum
  ## Psi_x is mult-dimensional, x can be of two dimesnion.

  Psi_x_minus <- function(x1, mu,  FIM,  x, w,  answering){
    ## mu and answering are only to avoid having another if when we want to check the maximum of sensitivity function
    Out <- arg$Psi_x(x1 = x1, mu =  mu, FIM = FIM,  x = x, w = w, answering = answering)
    return(-Out)
  }




  if(length(arg$lx) == 1)
    Psi_x_plot <-  arg$Psi_x ## for PlotPsi_x

  # it is necessary to distniguish between Psi_x for plotiing and finding ELB becasue in plotting for models with two
  # explanatory variables the function should be defined as a function of x, y (x, y here are the ploints to be plotted)

  if(length(arg$lx) == 2)
    Psi_x_plot <- arg$Psi_xy

  #when length(lx) == 1, then Psi_x_plot = Psi_x
  ################################################################################



  ###############################################################################
  # required for finding the answering set for verification
  #if (length(lp) <= 2)
  optim_starting <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  n_independent){
    out <- optim2(fn = fn, lower = lower, upper = upper,
                  n.seg = control$n.seg,
                  q = c(x, w),
                  fixedpar = fixedpar, fixedpar_id = fixedpar_id,
                  n_independent= n_independent)
    minima <- out$minima
    counts <- out$counts
    return(list(minima =minima, counts = counts))
  }

  ############################################################################################################

  ############################################################################################################
  # defining the nloptr directL
  if (type != "optim_on_average")
    vertices_inner <- make_vertices(lower = lp, upper = up)
  optim_nloptr <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  n_independent, inner_maxeval = control$inner_maxeval){

    out_directL <- nloptr::directL(fn = fn,
                                   lower = lower,
                                   upper = upper,
                                   q = c(x, w),
                                   nl.info = FALSE,
                                   fixedpar = fixedpar, fixedpar_id = fixedpar_id,
                                   n_independent = n_independent,
                                   control=list(xtol_rel = .Machine$double.eps,
                                                maxeval = inner_maxeval))
    out <- find_on_points(fn = fn,
                          points = vertices_inner,
                          q = c(x, w),
                          fixedpar = fixedpar,
                          fixedpar_id = fixedpar_id,
                          n_independent = n_independent)

    counts <- out$count + out_directL$iter
    minima <- out$minima
    minima_directL <- c(out_directL$par, out_directL$value)
    minima <- rbind(minima, minima_directL)
    return(list(minima =minima, counts = counts))
  }
  ############################################################################################################

  ############################################################################################################
  # for locally optimal design we dont need an inner optimization
  optim_locally <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  n_independent){

    log_det <- fn(param = param_locally, q=c(x, w),  fixedpar = fixedpar,
                  fixedpar_id = fixedpar_id, n_independent = n_independent)
    counts <- 1
    minima <- matrix(c(param_locally, log_det), nrow = 1)
    return(list(minima =minima, counts = counts))
  }
  ############################################################################################################


  ############################################################################################################
  # for locally optimal design we dont need an inner optimization
  value_optim_on_average <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  n_independent){

    log_det <- fn(param = param_set, q=c(x, w),  fixedpar = fixedpar,
                  fixedpar_id = fixedpar_id, n_independent = n_independent)
    counts <- 1
    minima <- cbind(param_set, log_det)
    return(list(minima =minima, counts = counts))
  }
  ############################################################################################################


  ############################################################################################################
  ## set the optim_func that will be used for calculating the cost values
  if (type != "locally" && type != "optim_on_average"){
    if (control$inner_space == "continuous"){
      optim_func <- optim_nloptr
    }
    if (control$inner_space == "vertices"){
      vertices <- make_vertices(lower = lp, upper = up)
      optim_func <- function(fn, lower, upper, w, x,
                             fixedpar, fixedpar_id, n_independent){
        out <- find_on_points(fn = fn,
                              points = vertices,
                              q = c(x, w),
                              fixedpar = fixedpar,
                              fixedpar_id = fixedpar_id,
                              n_independent = n_independent)
        minima <- out$minima
        counts <- out$counts
        return(list(minima =minima, counts = counts))
      }
    }

    if (control$inner_space == "discrete"){
      #discrete_set <- control$param_set

      optim_func <- function(fn, lower, upper, w, x,
                             fixedpar, fixedpar_id, n_independent){
        out <- find_on_points(fn = fn,
                              points = discrete_set,
                              q = c(x, w),
                              fixedpar = fixedpar,
                              fixedpar_id = fixedpar_id,
                              n_independent = n_independent)

        minima <- out$minima
        counts <- out$counts
        return(list(minima =minima, counts = counts))
      }
    }

  }
  if (type == "locally")
    optim_func <- optim_locally
  if(type == "optim_on_average")
    optim_func <- value_optim_on_average
  ############################################################################################################

  ############################################################################################################
  ## x_id, w_id are the index of x and w in positions
  #cost_id is the index of
  ## in symmetric case the length of x_id can be one less than the w_id if the number of design points be odd!
  if (control$sym)
    x_id <- 1:floor(arg$k/2) else
      x_id <- 1:(arg$k * n_independent)
  if (!control$equal_weight)
    w_id <- (x_id[length(x_id)] + 1):length(arg$ld) else
      w_id <- NA

  ###column index of cost in  matrix output of the inner problem
  if (type != "optim_on_average")
    CostColumnId <- length(lp) + 1 else
      CostColumnId <- dim(arg$param)[2] + 1


  ## warning: not the lp withot fixed param

  ######################################################################################################
  ## whenever Calculate_Cost is used, the fixed_arg list should be passed to
  ## fixed argumnet for function Calculate_Cost
  fixed_arg = list(x_id = x_id,
                   w_id = w_id,
                   sym = control$sym ,
                   sym_point = control$sym_point,
                   CostColumnId = CostColumnId,
                   crfunc2 = crfunc2,
                   lp = lp, ## NULL gets here
                   up = up, ## NULL gets here for optim_on_average
                   fixedpar = fixedpar,
                   fixedpar_id = fixedpar_id,
                   optim_func = optim_func,
                   n_independent = n_independent,
                   type = type,
                   equal_weight = control$equal_weight,
                   k = arg$k)
  if (type == "optim_on_average")
    fixed_arg$param <- arg$param

  ########################################################################################

  #################################################################################################
  # Initialization when evol is NULL
  #################################################################################################
  if (is.null(evol)){
    ## set the old seed if call is from minimax
    if (!is.null(control$rseed))
      set.seed(control$rseed)

    msg <- NULL
    revol_rate <- control$revol_rate
    maxiter <- iter
    totaliter <- 0
    #evol <- list()
    min_cost <- c() ## cost of the best imperialists
    mean_cost <- c() ## mean cost of all imperialists
    check_counter <- 0 ## counter to count the check
    total_nlocal  <- 0 ## total number of successful local search
    if (!control$lsearch)
      total_nlocal <- NA
    total_nrevol <- 0 ## total number of successful revolution
    total_nimprove <- 0 ##total number of improvements due to assimilation

    ############################################## Initialization for ICA
    InitialCountries <- GenerateNewCountry(NumOfCountries = control$ncount,
                                           lower = arg$ld,
                                           upper = arg$ud,
                                           sym = control$sym,
                                           w_id = w_id,
                                           x_id = x_id,
                                           n_independent= n_independent,
                                           equal_weight = control$equal_weight)


    if (!is.null(arg$initial))
      InitialCountries[1:dim(arg$initial)[1], ] <- arg$initial
    InitialCost <- vector("double", control$ncount)
    temp <- Calculate_Cost (mat = InitialCountries, fixed_arg = fixed_arg)
    total_nfeval <-  temp$nfeval
    InitialCost <-  temp$cost
    inparam <- temp$inner_optima
    ## waring inparam for optim_on_average does not have any meaning!
    temp <- NA # safety
    ##Now we should sort the initial countries with respect to their initial cost
    SortInd <- order(InitialCost)
    InitialCost <- InitialCost[SortInd] # Sort the cost in assending order. The best countries will be in higher places
    InitialCountries <- InitialCountries[SortInd,, drop = FALSE] #  Sort the population with respect to their cost. The best country is in the first column
    inparam <- inparam[SortInd, , drop = FALSE]
    # creating empires
    Empires <- CreateInitialEmpires(sorted_Countries = InitialCountries,
                                    sorted_Cost = InitialCost,
                                    Zeta = control$zeta,
                                    NumOfInitialImperialists = control$nimp,
                                    NumOfAllColonies = (control$ncount -control$nimp),
                                    sorted_InnerParam = inparam)
    best_imp_id<- 1 ## the index of list in which best imperialists is in.

    ########################################################################
  }
  #################################################################################################

  #################################################################################################
  # when we are updating the object for more number of iterations
  #################################################################################################
  if (!is.null(evol)){
    ## reset the seed!
    if (exists(".Random.seed")){
      GlobalSeed <- get(".Random.seed", envir = .GlobalEnv)
      #if you call directly from iterate and not minimax!
      on.exit(assign(".Random.seed", GlobalSeed, envir = .GlobalEnv))
    }

    msg <- object$best$msg
    prev_iter <- length(evol) ##previous number of iterationst
    maxiter <- iter + prev_iter
    totaliter <- prev_iter
    mean_cost <- sapply(1:(totaliter), FUN = function(j) evol[[j]]$mean_cost)
    min_cost <- sapply(1:(totaliter), FUN = function(j) evol[[j]]$min_cost)
    Empires <- object$empires


    check_counter <- arg$updating$check_counter
    total_nfeval <- object$alg$nfeval
    total_nlocal <-  object$alg$nlocal
    total_nrevol <-object$alg$nrevol
    total_nimprove <-object$alg$nimprove

    imp_cost <- round(sapply(object$empires, "[[", "ImperialistCost"), 12)
    best_imp_id<- which.min(imp_cost)
    revol_rate <-  arg$updating$revol_rate
    ##updating the random seed
    if (!is.null(control$rseed)){
      do.call("RNGkind",as.list(arg$updating$oldRNGkind))  ## must be first!
      assign(".Random.seed", arg$updating$oldseed , .GlobalEnv)
    }
  }
  ##########################################################################
  space_size <- arg$ud - arg$ld
  continue = TRUE

  #################################################################################################################
  ### start of the while loop until continue == TRUE
  #################################################################################################################
  while (continue == TRUE){
    totaliter <- totaliter + 1
    check_counter <- check_counter + 1
    # if (totaliter == 1058)
    #   browser()
    revol_rate <- control$damp * revol_rate
    ## revolution rate is increased by damp ration in every iter

    ###############################################################################################
    ################################################################# for loop over all empires[ii]
    for(ii in 1:length(Empires)){

      ########################################## local search is only for point!
      if (control$lsearch){
        LocalSearch_res <- LocalSearch (TheEmpire =  Empires[[ii]],
                                        lower = arg$ld,
                                        upper = arg$ud,
                                        l = control$l,
                                        fixed_arg = fixed_arg)

        Empires[[ii]] <- LocalSearch_res$TheEmpire
        total_nfeval <- total_nfeval + LocalSearch_res$nfeval
        total_nlocal <- total_nlocal + LocalSearch_res$n_success
      }
      ##########################################################################

      ############################################################## Assimilation
      temp5 <- AssimilateColonies2(TheEmpire = Empires[[ii]],
                                   AssimilationCoefficient = control$assim_coeff,
                                   VarMin = arg$ld,
                                   VarMax = arg$ud,
                                   ExceedStrategy = "perturbed",
                                   sym = control$sym,
                                   AsssimilationStrategy = control$assim_strategy,
                                   MoveOnlyWhenImprove = control$only_improve,
                                   fixed_arg = fixed_arg,
                                   w_id = w_id,
                                   equal_weight = control$equal_weight)
      ##Warning: in this function the colonies position are changed but the imperialist and the
      ##cost functions of colonies are not updated yet!
      ##they will be updated after revolution
      Empires[[ii]] <- temp5$TheEmpire
      total_nfeval <- total_nfeval + temp5$nfeval
      total_nimprove <-  total_nimprove + temp5$nimprove
      ##########################################################################

      ############################################################### Revolution
      temp4 <- RevolveColonies(TheEmpire = Empires[[ii]],
                               RevolutionRate = revol_rate,
                               NumOfCountries = control$ncount,
                               lower = arg$ld,
                               upper = arg$ud,
                               sym = control$sym,
                               sym_point = control$sym_point,
                               fixed_arg = fixed_arg,
                               w_id = w_id,
                               equal_weight = control$equal_weight)
      Empires[[ii]] <- temp4$TheEmpire
      total_nrevol <- total_nrevol + temp4$nrevol
      total_nfeval <- total_nfeval + temp4$nfeval
      ############################################################
      Empires[[ii]] <- PossesEmpire(TheEmpire = Empires[[ii]])

      ##after updating the empire the total cost should be updated
      ## Computation of Total Cost for Empires
      Empires[[ii]]$TotalCost <- Empires[[ii]]$ImperialistCost + control$zeta * mean(Empires[[ii]]$ColoniesCost)

    }
    ############################################################ end of the loop for empires [[ii]]
    ###############################################################################################

    #################################################### Uniting Similiar Empires
    if (length(Empires)>1){

      Empires <- UniteSimilarEmpires(Empires = Empires,
                                     Zeta = control$zeta,
                                     UnitingThreshold = control$uniting_threshold,
                                     SearchSpaceSize = space_size)
    }
    ############################################################################
    # zeta is necessary to update the total cost!
    Empires <- ImperialisticCompetition(Empires = Empires, Zeta = control$zeta)

    ############################################################## save the seed
    # we get the seed here because we dont know if in cheking it wil be chaned
    #we save the seed when we exit the algorithm
    oldseed <- get(".Random.seed", envir = .GlobalEnv)
    oldRNGkind <- RNGkind()
    ############################################################################

    ############################################################################
    # extracing the best emperor and its position
    imp_cost <- round(sapply(Empires, "[[", "ImperialistCost"), 12)
    min_cost[totaliter] <- switch(type, "minimax" = min(imp_cost), "standardized" = -min(imp_cost),
                                  "locally" =  min(imp_cost), "optim_on_average" = min(imp_cost))
    mean_cost[totaliter] <- switch(type, "minimax" = mean(imp_cost), "standardized" = -mean(imp_cost),
                                   "locally" = mean(imp_cost), "optim_on_average" = mean(imp_cost))

    best_imp_id <- which.min(imp_cost) ## which list contain the best imp
    if (!control$equal_weight)
      w <- Empires[[best_imp_id]]$ImperialistPosition[, w_id] else
        w <- w_equal
    x <- Empires[[best_imp_id]]$ImperialistPosition[, x_id]
    inparam <- Empires[[best_imp_id]]$ImperialistInnerParam
    if (length(lp)==1)
      inparam <- t(inparam)


    ##modifying the answering set if there is any fixed parameters.
    ## does not applicable for locally and optim_on_average
    if (any(!is.na(fixedpar))){
      fix_inparam <- c(fixedpar, inparam)
      NumOfParam <- 1:length(fix_inparam)
      inparam <- fix_inparam[order( c(fixedpar_id, setdiff(NumOfParam, fixedpar_id)))]
    }

    if (control$sym){
      x_w <- ICA_extract_x_w(x = x, w = w, sym_point = control$sym_point)
      x <- x_w$x
      w <- x_w$w
    }
    ##sort Point
    if (n_independent == 1){
      w <- w[order(x)]
      x <- sort(x)
    }
    ############################################################################

    ################################################################ print trace
    if (control$trace){
      if (type != "locally" && type != "optim_on_average")
        cat("\nICA iter:", totaliter, "\npoints:", x, "\nweights: ", w, "\nparam: ",
            inparam, "\nbest criterion value: ", min_cost[totaliter],"\n") else
              cat("\nICA iter:", totaliter, "\npoints:", x, "\nweights: ", w,
                  "\nbest criterion value: ", min_cost[totaliter],"\n")
    }
    ############################################################################

    if ( min_cost[totaliter] == 1e-24)
      warning("Computational issue! maybe the design is singular!\n")

    ################################################################### continue
    if (totaliter ==  maxiter){
      continue <- FALSE
      convergence = "maxiter"
    }
    if(length(Empires) == 1 && control$stop_rule == "one_empire"){
      continue <- FALSE
      convergence = "one_empire"
    }
    ## the continue also can be changed in check
    ############################################################################

    ################################################################################# plot_cost
    if (control$plot_cost) {
      #ylim for efficiency depends on the criterion type because
      ylim = switch(type,
                    "minimax" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2),
                    "standardized" = c(min(mean_cost[1:(totaliter)]) - .07, max( min_cost) + .2),
                    "locally" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2))
      PlotEffCost(from = 1,
                  to = (totaliter),
                  AllCost = min_cost, ##all criterion up to now (all cost function)
                  UserCost = NULL,
                  DesignType = type,
                  col = line_col[1],
                  xlab = "Iteration",
                  ylim = ylim,
                  lty = 1,
                  title1 = title1,
                  plot_main = TRUE)
      ICAMean_line <-  mean_cost[1:(totaliter)]
      lines(x = 1:(totaliter),
            y = ICAMean_line,
            col = line_col[2], type = "s", lty = 5)
      legend(x = legend_place,  legend = legend_text,lty = c(1,5, 3), col = line_col, xpd = TRUE, bty = "n")
    }
    ############################################################################

    ####################################################################################
    #check the equivalence theorem and find ELB
    ####################################################################################
    ## we check the quvalence theorem in the last iteration anyway. but we may not plot it.
    if (check_counter == control$equivalence_every || (check_last && !continue)){
      check_counter <- 0
      if (type != "locally" && type != "optim_on_average"){
        ########################################################################
        # find all local minima
        output <- optim_starting(fn= crfunc2,
                                 lower=lp,
                                 upper = up,
                                 w = w,
                                 x = x ,
                                 fixedpar = fixedpar,
                                 fixedpar_id = fixedpar_id,
                                 n_independent = n_independent)

        # total_nfeval <- total_nfeval + output$counts
        all_optima <- output$minima
        ########################################################################

        ########################################################################
        ### find the nearest elements in all_optima to construct the answering set larter
        near_ind <- find_nearest(values = all_optima[, dim(all_optima)[2]],
                                 tol = control$answering_merg_tol,
                                 compare_with_minimum = TRUE)
        ########################################################################

        ########################################################################
        ## add the columns of fixedpar to all local minima
        if (any(!is.na(fixedpar))){
          ## warnings: all_optima also contain the cost values in last column!
          ## create the columns of fixedpar
          fixedpar_col <- matrix(fixedpar, dim(all_optima)[1], length(fixedpar), byrow = TRUE)
          ## how many columns the all_optima with the fixed param will have
          all_dim <- 1:(dim(fixedpar_col)[2] + dim(all_optima)[2])
          all_optima <- cbind(all_optima, fixedpar_col)[, order(c(setdiff(all_dim, fixedpar_id), fixedpar_id)), drop = FALSE]
          #Modifying the all local_minima
        }
        ########################################################################

        # construct the answering set and alos answering set cost
        answering <- all_optima[near_ind, , drop = FALSE]
        answering_cost <- answering[, dim(answering)[2], drop = FALSE]
        if (type == "minimax")
          answering_cost <- -answering_cost
        answering <- answering[, -dim(answering)[2], drop = FALSE]
        ########################################################################

        ########################################################################
        ## we save all minima as well
        all_optima_cost <- all_optima[, dim(all_optima)[2], drop = FALSE]
        if (type == "minimax")
          all_optima_cost <- -all_optima_cost ## because we multiplied the inner problme by minus for minimax optimal design for directL
        all_optima <- all_optima[, -dim(all_optima)[2], drop = FALSE]
        ########################################################################

        ########################################################################
        # find the measure
        mu <- find_measure (npar = npar, x = x, w = w,
                            answering = answering,
                            FIM = arg$FIM,
                            Psi_Mu = arg$Psi_Mu)$mu
        ########################################################################
      }else{
        ## if type is locally
        if (type == "locally"){
          answering <- matrix(param_locally, nrow = 1) ## we need it for find measure and check. they use answering set
          mu <- 1 # we need it for check and plot
        }
        if (type == "optim_on_average"){
          answering <- arg$param
          mu <- arg$prior
        }


        answering_cost <- all_optima <- all_optima_cost <- NA
      }
      ##########################################################################
      # find the maximum of sensitivityfunction
      OptimalityCheck <- directL(fn = Psi_x_minus,
                                 lower = arg$lx,
                                 upper = arg$ux,
                                 mu = mu,
                                 answering = answering,
                                 x = x,
                                 w = w,
                                 FIM = arg$FIM,
                                 nl.info = FALSE,
                                 control=list(xtol_rel=.Machine$double.eps,
                                              maxeval = control$maxeval_equivalence))


      ##sometimes the optimization can not detect maximum  in the bound, so here we add the cost values on the bound
      vertices_outer <- make_vertices(lower = arg$lx, upper = arg$ux)
      check_vertices <- find_on_points(fn = arg$Psi_x,
                                       points = vertices_outer,
                                       mu = mu,
                                       answering = answering,
                                       x = x,
                                       w = w,
                                       FIM = arg$FIM)
      vertices_val <- check_vertices$minima[, dim(check_vertices$minima)[2]]
      ## minus because for optimality check we minimize the minus sensitivity function
      max_deriv <- c(-OptimalityCheck$value, vertices_val)
      max_deriv <- max(max_deriv)
      ##########################################################################

      # Efficiency lower bound
      ELB <- npar/(npar + max_deriv)

      ###Point_mat <- matrix(x, ncol = n_independent, nrow = arg$k)

      GE_confirmation <- (ELB >= control$stoptol)
      ##########################################################################
      # print trace that is related to checking
      if (control$trace)
        cat("maximum of sensitivity:", max_deriv, "\nefficiency lower bound (ELB):", ELB, "\n")
      ##########################################################################

      #if (n_independent == 1){
      if (GE_confirmation && control$stop_rule == "equivalence"){
        continue <- FALSE
        convergence <- "equivalence"
      }
      #}
      ########################################################### plot derivative
      ## plot sensitivityfunction
      if (control$plot_deriv)
        PlotPsi_x(lower =    arg$lx ,
                  upper =      arg$ux ,
                  Psi_x = Psi_x_plot,
                  FIM  = arg$FIM,
                  mu = mu,
                  x = x,
                  w = w,
                  answering = answering)
      ##########################################################################
    }else
      max_deriv <- answering <- answering_cost <-all_optima <- all_optima_cost  <- mu <- ELB <- NA
    ####################################################################### end of check
    #  if (check_counter == control$equivalence_every || (check_last && !continue)) ##########
    ####################################################################################

    if (type == "locally" || type == "optim_on_average"){
      answering <- NA # now we dont need answering. We required it before for checking so we set it to NA
      mu <- 1
    }

    ####################################################################### save
    evol[[totaliter]] <- list(iter = totaliter,
                              x = x,
                              w = w,
                              min_cost = min_cost[totaliter],
                              mean_cost = mean_cost[totaliter],
                              all_optima = all_optima,
                              all_optima_cost = all_optima_cost,
                              answering = answering,
                              answering_cost = answering_cost,
                              mu = mu,
                              max_deriv = max_deriv,
                              ELB = ELB)

    if (type != "locally" && type != "optim_on_average"){
      evol[[totaliter]]$param = inparam
    } else
      evol[[totaliter]]$param = NA
    ############################################################################

    ################################################################ print trace
    if (control$trace){
      cat("total local search:", total_nlocal, "\n")
      cat("total revolution:", total_nrevol, "\n")
      if (control$only_improve)
        cat("total improve:", total_nimprove, "\n")
    }
    ############################################################################
  }
  #################################################################################################################
  ### end of the while loop over continue == TRUE
  #################################################################################################################

  if (!control$only_improve)
    total_nimprove <- NA

  #if (ELB >= control$stoptol && control$stop_rule == "equivalence")
  #  convergence = "equivalence" else

  ##############################################################################
  # check the appropriateness of the maxeval
  # if (type != "locally" & type != "optim_on_average" & control$inner_space == "continuous"){
  #   if (control$check_inner_maxeval){
  #     check_temp <- check_maxeval(fn = crfunc2, lower = lp, upper = up, maxeval = control$inner_maxeval,
  #                                 fixedpar = fixedpar, fixedpar_id = fixedpar_id, n_independent = n_independent, q = c(x, w))
  #     msg <- check_temp$msg
  #   }else
  #     msg <- NULL
  # }

  ##################
  msg <- NULL
  ##############################################################################

  ######################################################################## saving
  ## we add the following to arg becasue dont want to document it in Rd files
  # updating parameters
  object$arg$updating$check_counter <- check_counter
  object$arg$updating$oldseed <- get(".Random.seed", envir = .GlobalEnv)
  object$arg$updating$oldRNGkind <-  RNGkind
  object$arg$updating$revol_rate = revol_rate ## different from revolrate

  object$evol <- evol
  object$empires <- Empires
  # object$best <- list(
  #   #empires = Empires,
  #   iter = totaliter,
  #   x = evol[[totaliter]]$x,
  #   w = evol[[totaliter]]$w,
  #   cr = evol[[totaliter]]$min_cost,
  #   all_optima = evol[[totaliter]]$all_optima,
  #   all_optima_cost = evol[[totaliter]]$all_optima_cost,
  #   answering = evol[[totaliter]]$answering,
  #   answering_cost = evol[[totaliter]]$answering_cost,
  #   max_deriv = evol[[totaliter]]$max_deriv,
  #   ELB = evol[[totaliter]]$ELB,
  #   mu = evol[[totaliter]]$mu,
  #   nfeval = total_nfeval,
  #   nlocal = total_nlocal,
  #   nrevol = total_nrevol,
  #   nimprove = total_nimprove,
  #   convergence = convergence,
  #   msg = msg)


  object$alg <- list(
    nfeval = total_nfeval,
    nlocal = total_nlocal,
    nrevol = total_nrevol,
    nimprove = total_nimprove,
    convergence = convergence)
    #msg = msg

  ## so object is 'res', 'arg' and 'evol'
  ## arg has a list named update as well
  ###### end of saving
  ##############################################################################
  return(object)

}


