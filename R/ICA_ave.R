#' Imperialist Competitive Algorithm to find optimum on-the-average designs based on the least favorable distribution.
#'
#' For nonlinear models, if a set of possible values of the vector parameters \eqn{\bold{\theta}} is available, the optimum-on-the-average designs
#' can be applied. In this approach the criterion is evaluated at plausible parameter values and weighted by a probability measure
#' \eqn{\pi}, the measure having support in the parameter space \eqn{\Theta}. The resulted weighted criterion is thus
#' \deqn{B(\xi, \pi) = \int_{\Theta}\Psi(\xi, \theta)d(\pi(\theta)).}{B(\xi, \Pi) = \int_{\Theta}\Psi(\xi, \theta)d(\pi(\theta)).}
#' For current version \eqn{\Psi(\xi, \theta)} is equal to \eqn{|M(\xi, \theta)|}.
#'
#'
#'
#' @param fimfunc  the name of the FIM from available FIM functions in the package as a character string or the user-written function that returns the FIM as a \code{matrix}. See "Details" in \link{mica}.
#' @param lx lower bound of the design space \eqn{\chi}
#' @param ux upper bound of the design space \eqn{\chi}
#' @param prior a vector of the probability measure \eqn{\pi}.
#' @param param a matrix for set of parameters, i.e. support of \eqn{\pi}. Every row is is a vector of values of a parameter.
#' The number of its rows must be equal to the length of \code{prior}.
#' @param iter maximum number of iterations
#' @param k number of design (support) points. Must be larger than the number of model parameters \eqn{p} to avoid singularity of the FIM.
#' @param control a list of contains the tuning parameters of ICA and the design problem. See "Details" of \code{\link{mica}}.
#' @param initial a matrix of the  initial designs that will be used as initial countries in ICA.
#'  Every row is a design and concatenation of \code{x} and \code{w}.  Default is \code{NULL}. See "Details" of \code{\link{mica}}.
#' @param ... further arguments to be passed to the FIM function given by \code{fimfunc}.
#'
#' @return
#' an object of class "ICA". See "Value" in \code{\link{mica}}.
#' @examples
#' \dontrun{
#' test <- ave(fimfunc = "FIM_logistic",
#'                        lx = -5, ux = 5, prior = rep(1/4, 4),
#'                        param = matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2),
#'                        iter = 200, k = 3)
#'
#' plot(test)
#' print(test)
#'################################################################################
#'## using equivalence theorem as stopping rule. Can be applied in H-algorithm
#' test <- ave (fimfunc = "FIM_logistic",
#'                          lx = -5, ux = 5, prior = rep(1/4, 4),
#'                          param = matrix(c(0.5, 1.5, 0.5, 1.5, 4.0, 4.0, 5.0, 5.0), 4, 2),
#'                          iter = 200, k =3,
#'                          control = list(stop_rule = "equivalence",
#'                          stoptol = .9995, equivalence_every = 100))
#'}
#'
#' @export

ave <- function(fimfunc,
                           lx,
                           ux,
                           prior,
                           param,
                           iter,
                           k,
                           control = list(),
                           initial = NULL,
                           ...) {

  #############################################################################
  #### common with mfw
  if (missing(fimfunc))
    stop("\"fimfunc\" is missing")
  if (missing(lx))
    stop("\"lx\" is missing")
  if (missing(ux))
    stop("\"ux\" is missing")
  if (missing(param))
    stop("\"param\" is missing")
  if (missing(prior))
    stop("\"prior\" is missing")
  if (missing(iter))
    stop("\"iter\" is missing")
  if (missing(k))
    stop("\"k\" is missing")
  if (!is.function(fimfunc) && !is.character(fimfunc))
    stop(" \"fimfunc\" can be either \"character\" or \"function\"")
  if (length(lx) != length(ux))
    stop("Length of \"lx\" is not equal to length of \"ux\"")
  if (length(prior) != dim(param)[1])
    stop("length of \"prior\" is not equal to the number of rows of \"param\"")
  if (!is.numeric(k) || (k %% 1) != 0 || k <= 0)
    stop("\"k\" must be a positive integer number")
  if (k < dim(param)[2])
    stop("\"k\" must be larger than the number of parameters to avoid singularity")
  if (!is.numeric(iter) || (iter %% 1) != 0 || iter <= 0)
    stop("\"iter\" must be a positive integer number")
  #############################################################################
  type = "D"
  if (missing(type))
    stop("'type' is missing")
  if (!(type %in% c("D", "standardized_D")))
    stop ("'type' must be \"D\" or \"standardized_D\"")

  #############################################################################
  ## ICA tuning parameters
  if (is.null(control$ncount))
    control$ncount <- 40
  if (is.null(control$nimp))
    control$nimp <- control$ncount / 10 ## 10 percent of the countries are imperialists
  if (control$ncount - control$nimp <= control$nimp)
    stop(
      "Number of colonies is less than the number of imperialists. Please increase the number of countries or decrease the number of imperialists."
    )
  if (is.null(control$assim_coeff))
    control$assim_coeff <- 4
  if (is.null(control$revol_rate))
    control$revol_rate <- .3
  if (is.null(control$damp))
    control$damp <-
      .99 ## the revolution rate decreases by this rate
  if (is.null(control$zeta))
    control$zeta <-
      .1 ### the role of colonist in calculating the total cost of empire. page 1224 optimum design for skeletal structures by A. Kaveh and Talatahari
  if (is.null(control$uniting_threshold))
    control$uniting_threshold <- .02
  if (is.null(control$assim_strategy))
    control$assim_strategy <- "PICA"
  if (is.null(control$lsearch))
    control$lsearch <- 2
  if (is.null(control$l))
    control$l <- 2
  if (is.null(control$only_improve))
    control$only_improve <- TRUE
  #############################################################################

  #############################################################################
  # some mica functionality for stopping rules
  if (is.null(control$stop_rule))
    control$stop_rule <- "maxiter" # "equivalence" or "one_empire"
  if (is.null(control$stoptol))
    control$stoptol <- .99
  if (is.null(control$equivalence_every))
    control$equivalence_every <- 200
  #############################################################################

  #############################################################################
  ## about the design space and weights:
  if (is.null(control$equal_weight))
    control$equal_weight <- FALSE
  ## if TRUE, then ld and ud does not have the lower bound and upper bound for the weights.
  ## In this case, the countries are only the points and not weights
  if (is.null(control$sym))
    control$sym <- FALSE
  if (control$sym && is.null(control$sym_point))
    stop("the symetric point should be provided by 'sym_point' in 'control'")
  # for logistic mean(c(lp[1], up[1]))
  #sym_poitn will remain NULL
  if (length(lx) != 1 && control$sym)
    ## the symetrix is only for model with one indepenent variable
    stop("currently symetric property only can be applied to models with one variable")
  if (control$equal_weight && control$sym)
    stop("symmetric property does not work when only equal-weighted design is requested")
  #############################################################################

  # #############################################################################
  ## inner space
  if (is.null(control$inner_space))
    control$inner_space <- "none" # or "vertices"

  # if (!control$inner_space %in% c("continuous", "vertices", "discrete"))
  #   stop("inner space can be \"continuous\", \"vertices\" or \"discrete\"")
  # if (control$inner_space == "discrete" && is.null(control$param_set))
  #   stop("'param_set' must be given")
  # if (is.null(control$inner_maxeval))
  #   control$inner_maxeval <- 600
  # if (is.null(control$check_inner_maxeval))
  #   control$check_inner_maxeval <- TRUE
  # #############################################################################

  #############################################################################
  ## plots and trace and seed
  if (is.null(control$plot_cost)){
    control$plot_cost <- FALSE
  }
  if (is.null(control$plot_deriv))
    control$plot_deriv <- TRUE
  if (is.null(control$trace))
    control$trace <- TRUE
  if (is.null(control$rseed))
    control$rseed <- NULL
  #############################################################################

  #############################################################################
  ## for gosolnp
  if (is.null(control$gosolnp_trace))
    control$gosolnp_trace <- FALSE
  if (is.null(control$gosolnp_n.restarts))
    control$gosolnp_n.restarts <- 1
  if (is.null(control$gosolnp_n.sim))
    control$gosolnp_n.sim <- 800

  ## list
  control_gosolnp <- list(
    gosolnp_trace  = control$gosolnp_trace,
    gosolnp_n.restarts = control$gosolnp_n.restarts,
    gosolnp_n.sim = control$gosolnp_n.sim,
    gosolnp_rseed =  NULL
  )
  #############################################################################

  #############################################################################
  ## to check equivalence theorem
  if (is.null(control$n.seg))
    control$n.seg <- 4
  if (is.null(control$maxeval_equivalence))
    control$maxeval_equivalence <- 6000 ## maxiter for directL that is used to find the maximum of the derivative function
  # not by user
  control$answering_merg_tol <- .005
  #############################################################################

  #############################################################################
  ### do not change the seed
  if (exists(".Random.seed")) {
    OldSeed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", OldSeed, envir = .GlobalEnv))
  }
  #############################################################################

  #############################################################################
  # setting the Fisher information matrix
  #if character, set fimchar and find the appropriate FIM function.
  if (is.character(fimfunc)) {
    fimchar <- fimfunc
    fimfunc <- check_lp_fimchar(fimchar = fimchar, lp = param[1, ]) ## to only check the length of parameters
  }else
    fimchar <- "unknown"

  # to handle ...
  fimfunc2 <- function(x, w, param) {
    out <- fimfunc(x = x, w = w, param = param,...)
    return(out)
  }
  #############################################################################





  # global variables needed for the definition of crfunc
  n_independent <- length(lx)
  npar <- dim(param)[2]

  ########
  type <- paste("on_average", type, sep = "_")

  #############################################################################
  ## creating the criterion function
  ## prior and param_set are from global environment
  if (type == "on_average_D"){ ## D-criterion
    crfunc_on_average_D <- function(param, q, n_independent) {
      lq <- length(q) # q is the design points and weights
      pieces <- lq / (n_independent + 1)
      x_ind <- 1:(n_independent * pieces)
      w_ind <- (x_ind[length(x_ind)] + 1):lq
      x <- q[x_ind]
      w <- q[w_ind]
      on_average_crfunc <- sum(
        sapply(1:length(prior), FUN= function(j)
          prior[j] * -det2(fimfunc2(x = x, w = w, param = param[j, ]), logarithm = TRUE))
      ) + 5000 * (sum(w) - 1) ^ 2   ## -(-det+pen) = det-pen
      return(on_average_crfunc)
    }
  }

  if (type == "on_average_standardized_D"){

    crfunc_locally <- function(param, q, n_independent) {

      lq <- length(q)
      pieces <- lq / (n_independent + 1)
      x_ind <- 1:(n_independent * pieces)
      w_ind <- (x_ind[length(x_ind)] + 1):lq
      x <- q[x_ind]
      w <- q[w_ind]
      locally_det <-
        -det2(fimfunc2(x = x, w = w, param = param), logarithm = TRUE) + 5000 *
        (sum(w) - 1) ^ 2
      if (locally_det == -1e24)
        ## becuase for locally the 'locally_det' will be -1e24 and spoil the algorithm!
        locally_det <- 1e24
      return(locally_det)
    }

    # ld <- c(rep(lx, npar), rep(0, npar))
    # ud <- c(rep(ux, npar), rep(1, npar))
    # quiqer but without equivalence theorem
    # locally_res <- sapply(1:length(prior), function(j)Rsolnp::gosolnp(fun = crfunc_locally, LB = ld, UB = ud,
    #                                                    n.restarts = control_gosolnp$gosolnp_n.restarts,
    #                                                    n.sim = control_gosolnp$gosolnp_n.sim,
    #                                                    control= list(trace = control_gosolnp$gosolnp_trace),
    #                                                    param = param[j, ], n_independent = n_independent))
    n_prior <- length(prior)
    locally_res <- sapply(1:n_prior, function(j)mica(fimfunc = fimfunc2, lx = lx, ux = ux,
                                                     lp = param[j, ], up = param[j, ], iter = 2000, k = npar,
                                                     type = "locally",
                                                     control = list(plot_cost = FALSE, plot_deriv = FALSE, trace = FALSE,
                                                                    stop_rule = "equivalence", equal_weight = TRUE,
                                                                    stop_tol = .999,  equivalence_every = 100)))
    locally_val <- rep(NA,  n_prior)
    ELB_val <- rep(NA,  n_prior)
    locally_design <- matrix(NA, ncol = dim(param)[1], nrow = n_prior)

    for (i in 1:n_prior){
      locally_val[i] <- locally_res[, i]$evol[[length(locally_res[, i]$evol)]]$min_cost
      ELB_val[i] <- locally_res[, i]$evol[[length(locally_res[, i]$evol)]]$ELB
      locally_design[i, ] <- c(locally_res[, i]$evol[[length(locally_res[, i]$evol)]]$x,
                               locally_res[, i]$evol[[length(locally_res[, i]$evol)]]$w)
    }
    locally_val <- exp(-locally_val)
    crfunc_on_average_standardized_D <- function(param, q, n_independent) {
      lq <- length(q) # q is the design points and weights
      pieces <- lq / (n_independent + 1)
      x_ind <- 1:(n_independent * pieces)
      w_ind <- (x_ind[length(x_ind)] + 1):lq
      x <- q[x_ind]
      w <- q[w_ind]
      on_average_crfunc <- sum(
        sapply(1:length(prior), FUN= function(j)
          -prior[j] * det2(fimfunc2(x = x, w = w, param = param[j, ]), logarithm = FALSE)/locally_val[j])
      ) + 5000 * (sum(w) - 1) ^ 2   ## -(-det+pen) = det-pen
      return(on_average_crfunc)
    }
  }
    #############################################################################


  if (type == "on_average_D")
    crfunc <- crfunc_on_average_D

  if (type == "on_average_standardized_D")
    crfunc <- crfunc_on_average_standardized_D
  #############################################################################
  # create the upper bound and lower bound for and the dimension of decision variables
  #x_length can be
  if (control$sym) {
    ## if the number of design points be odd then the middle point of the design shouyld be the symmetric point!
    ## the number of design point can be one less then w for odd k
    if (k %% 2 == 0) {
      w_length <-  k / 2
      x_length <-  k / 2
    }else{
      x_length <-  floor(k / 2)
      w_length <- floor(k / 2) + 1
    }

  }else{
    x_length  <-  k * n_independent
    w_length <-  k
  }
  # so if sym == TRUE then we have two possiblity:
  # the number of design point is odd, then the x is one element less than w and in the
  # the missed poiint is the symetric point!
  # but if number of design be even then the number length of x is equal to length of w
  ld <- c()
  ud <- c()
  for (i in 1:n_independent) {
    ld <- c(ld, rep(lx[i], x_length / n_independent))
    ud <- c(ud, rep(ux[i], x_length / n_independent))
  }

  #now we should add the wights!
  if (!control$equal_weight) {
    ld = c(ld, rep(0, w_length))
    ud = c(ud, rep(1, w_length))
  }
  #############################################################################


  #############################################################################
  ## dealing with initial countries if set by user
  if (!is.null(initial)) {
    # convert the intial to matrix if it is a vector!
    if (!is.matrix(initial))
      initial <- t(as.matrix(initial))
    initial <- round(initial, 8) #because it may produce strange error when cheking the lower upper bound!!!
    # check if the length is true!
    if (!is.na(initial) && dim(initial)[2] != length(ld))
      stop("The number of columns of 'initial' does not match with length of countries and should be", length(ld))
    #  we use round to protect the function from strange behaviour
    initial_out <-  sapply(1:dim(initial)[1], function(j) any(round(initial[j,], 5) > round(ud, 5)) || any(round(initial[j,], 5) < round(ld, 5)))
    if (any(initial_out))
      stop("The initial vaule(s) in row ", paste(which(initial_out), collapse = ", "), " are (is) out of bound.")
  }
  #############################################################################


  #############################################################################
  ## making the arg list
  ## the variables that will be added to control not by user, but by mica
  arg <-
    list(
      lx = lx, ux = ux, k = k, ld = ld, ud = ud,
      FIM = fimfunc2, crfunc = crfunc,
      initial = initial, control = control, type =type,
      prior = prior, param = param
    )
  ## updating will be added to arg in iterate
  ### sensitivity function required ofr cheking the equivalence theorem


  #############################################################################
  ### sensitivity function required ofr cheking the equivalence theorem
  arg$Psi_x <- Psi_x
  arg$Psi_Mu <- Psi_Mu
  if (length(lx) == 2)
    arg$Psi_xy <- Psi_xy
  ##  Psi_x works for both one, two and three dimensional
  ## but Psi_xy is needed for plotting becasue the function should have two arguments
  #############################################################################

  ICA_object <- list(arg = arg,
                     evol = NULL)
  class(ICA_object) <- c("list", "ICA")

  out <- iterate.ICA(object = ICA_object, iter = iter)
  return(out)

}

