#' Imperialist Competitive Algorithm to find multiple-objective optimal designs for the 4-parameter logistic models.
#'
#' The 4-parameter Hill model is \eqn{Y_i = f(D_i, a, b, c, d) + \epsilon_i}{Y_i = f(Di, a, b, c, d) + \epsiloni}
#' where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsiloni ~ N(0, \sigma2)},
#'  \deqn{f(D_i, a,b, c, d) = c + \frac{(d-c)(\frac{D_i}{a})^b}{1+(\frac{D_i}{a})^b}}{
#'  f(Di, a,b, c, d) = c + (d-c)(Di/a)^b/(1 + (1 + Di)^b)} is the mean response at dose \eqn{D_i}{Di}, \eqn{a} is the ED50,
#'  \eqn{d} is the upper limit of response, \eqn{c} is the lower limit of response and \eqn{b} denotes the Hill constant that
#'  control the flexibility in the slope of the response curve.\cr
#' Let \eqn{\bold{\theta} = (\theta_1, \theta_2, \theta_3, \theta_4)}{\bold{\theta} = (\theta1, \theta2, \theta3, \theta4)}. The Hill model
#'  may be re-parameterized as
#'  \deqn{f(x, \bold{\theta}) = \frac{\theta_1}{1 + exp(\theta_2 x + \theta_3)} + \theta_4,}{
#'  f(x, \bold{\theta})= \theta1/(1 + exp(\theta2*x + \theta3)) + \theta4,}
#'  which is sometimes referred to as \strong{4-parameter logistic model}. This form is equivalent to
#'  Hill model with \eqn{\theta_1 = d - c}{\theta1 = d - c}, \eqn{\theta_2 = - b}{\theta2 = - b},
#'   \eqn{\theta_3 = b\log(a)}{\theta3 = blog(a)}, \eqn{\theta_4 = c}{\theta4 = c}, \eqn{\theta_1 > 0}{\theta1 > 0},
#'   \eqn{\theta_2 \neq 0}{\theta2 not equal to 0}, and \eqn{-\infty < ED50 < \infty}, where \eqn{x_i = log(D_i) \in [-M, M]}{xi = log(Di) belongs to [-M, M]}
#'   for some sufficiently large value of \eqn{M}. The user can easily transform \eqn{a}, \eqn{b}, \eqn{c}, \eqn{d} and \eqn{D}
#'    to \eqn{\theta_1}{\theta1}, \eqn{\theta_2}{\theta2}, \eqn{\theta_3}{\theta3}, \eqn{\theta_4}{\theta4}
#'     and \eqn{x} to find the multiple-objective optimal design for the 4-parameter logisitc model with this function and then
#'     transform the design for 4-parameter logistic model  to the 4-parameter Hill model. See "Examples". \cr
#'
#'
#'
#' @param lx lower bound of the design space \eqn{\chi}.
#' @param ux upper bound of the design space \eqn{\chi}.
#' @param param initial guess of parameters \eqn{\bold{\theta} = (\theta_1, \theta_2, \theta_3, \theta_4)}{\bold{\theta} = (\theta1, \theta2, \theta3, \theta4)}.
#' @param iter maximum number of iterations.
#' @param k number of design (support) points. Must be larger than the number of model parameters \eqn{4} to avoid singularity of the FIM.
#' @param lambda user select weights, where \eqn{\lambda_1}{\lambda1} is the weight for estimating parameters,
#' \eqn{\lambda_2}{\lambda2} is the weignt for estimating median effective dose level (ED50), and \eqn{\lambda_3}{\lambda3} is the weight for estimating minimum effective dose level (MED).
#' @param delta numeric, predetermined clinically significant effect to define the MED.
#' When the dose-response relationship (slop \eqn{\theta_2}{\theta2}) is decreasing (increasing), the value of \code{delta} is negative (positive).
#' @param control a list of control parameters. See "Details" of \code{\link{mica}}.
#' @param initial initial a matrix of user intial countries or
#' a vector of a country that will be inserted  into the initial countries of ICA. See "Details" of \code{\link{mica}}.
#' @references
#' Hyun, S. W., & Wong, W. K. (2015). Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. The international journal of biostatistics, 11(2), 253-271.
#' @details
#' When \eqn{\lambda_1 > 0}{\lambda1 > 0}, then the number of support points \code{k} must at least be 4 to avoid singularity of the information matrix.\cr
#' This function is not suitable for finding the c-optimal designs for estimating 'MED' or 'ED50' and the results may not be stable.
#'  The reason is that for c-optimal the generalized inverese of Fisher information matrix unstable and depends
#'  on the tolerenace (probably because the matrix is exactly singular). \cr
#'
#'  The tolerance for finding the general inverse of Fisher information matrix is set to \code{.Machine$double.xmin}.
#'
#' @return
#' an object of class "ICA". See "Value" in \code{\link{mica}}.
#' @export
#' @examples
#'
#'## An example on how to create the design in Hyun and Wong (2015)
#'## An initial guess from Table 1:
#'
#'Theta1 <- c(1.563, 1.790, 8.442, 0.137)
#'
#'#########################################################
#'## Table 2 first row
#'# creating optimal design for estimating parameters: xi_D
#'res <- multica_4pl(lx = log(.001),
#'                    ux = log(1000),
#'                    param = Theta1,
#'                    k = 4,
#'                    lambda = c(1, 0, 0),
#'                    delta = -1,
#'                    iter = 150,
#'                    control = list(seed = 1366, plot_cost = TRUE))
#'\dontrun{
#'#######################################################
#'## finding multiple objective optimal design: example 1, Table 3
#'res1 <- multica_4pl(lx = log(.001),
#'                    ux = log(1000),
#'                    param = Theta1,
#'                    k = 4, lambda = c(0.05, 0.05, .90),
#'                    delta = -1, iter = 400,
#'                    control = list(seed = 1366))
#'
#'plot(res1)
#'
#'#######################################################
#'## finding multiple objective optimal design: example 2, Table 3
#'res2 <- multica_4pl(lx = log(.001),
#'                    ux = log(1000),
#'                    param = c(16.8, -1, 4.248, 22),
#'                    k = 4, lambda = c(0, 0.1, .9),
#'                    delta = 5, iter = 200,
#'                    control = list(seed = 1366))
#'plot(res2)
#'
#'##########################################################
#'## how to transfer from Hill model to 4-parameter logistic model
#'## parameters for Hill model
#'a <- 0.008949  # ED50
#'b <- -1.79 # Hill constant
#'c <- 0.137 # lower limit
#'d <- 1.7 # upper limit
#'D <- c(.001, 1000) ## dose in mg
#'## Hill_para is c(a, b, c, d)
#'res2 <- multica_4pl(lx = log(.001),
#'                    ux = log(1000),
#'                    param =  c(d - c, -b, b * log(a), c),
#'                    k = 4, lambda = c(0.05, 0.05, .90),
#'                    delta = -1, iter = 400,
#'                    control = list(seed = 1366))
#'exp(res2$evol[[length(res2$evol)]]$x) # dose level in mg
#'}
multica_4pl <- function(lx,
                        ux,
                        param,
                        iter,
                        k,
                        delta,
                        lambda,
                        control = list(),
                        initial = NULL) {


  ################################
  ### fimfunc
  fimfunc <- "FIM_logistic_4par"
  lp <- up <- param
  type <- "multiple_locally"
  multiple <- list(delta = delta, lambda = lambda)
  ## delta will be used when creating the multi objective senisitivity function and criterion
  #################################
  if (missing(fimfunc))
    stop("\"fimfunc\" is missing")
  if (missing(lx))
    stop("\"lx\" is missing")
  if (missing(ux))
    stop("\"ux\" is missing")
  if (missing(lp))
    stop("\"lp\" is missing")
  if (missing(up))
    stop("\"up\" is missing")
  if (missing(param))
    stop("\"param\" is missing")
  if (missing(iter))
    stop("\"iter\" is missing")
  if (missing(k))
    stop("\"k\" is missing")

  #################################################################################
  if (!is.function(fimfunc) && !is.character(fimfunc))
    stop(" \"fimfunc\" can be either \"character\" or \"function\"")
  #################################################################################

  # if (!is.character(fimfunc))
  #   stop(" \"fimfunc\" must be a character string")
  if (length(lx) != length(ux))
    stop("Length of \"lx\" is not equal to length of \"ux\"")
  if (length(lp) != length(up))
    stop("length of \"lp\" is not equal to length of \"up\"")
  if (!is.numeric(k) || (k %% 1) != 0 || k <= 0)
    stop("\"k\" must be a positive integer number")

  if (!type %in% c("minimax", "standardized", "locally", "multiple_locally"))
    stop("\"type\" must be \"minimax\", \"standardized\",  \"locally\",  \"multiple_locally\"")
  if (type == "locally" & !all(lp == up) || type == "multiple_locally" & !all(lp == up))
    stop("\"lp\" must be equal to \"up\" for locally (multiple-objective) optimal designs")
  if (!is.numeric(iter) || (iter %% 1) != 0 || iter <= 0)
    stop("\"iter\" must be a positive integer number")


  ############################################################### multiobjective optimal design checking

  if (type == "multiple_locally"){
    if (fimfunc == "FIM_logistic_4par"){

      if (is.null(delta) | is.null(lambda))
        stop("'delta' and 'lambda' must be given")
      if (param[2]>0 && delta > 0 || param[2]<0 && delta < 0)
        stop("'delta' must be negative when theta2 > 0, or positive when theta2 < 0")
      if (lambda[1] >0 & k < length(lp))
        stop("\"k\" must be larger than the number of parameters (4) to avoid exact singularity")
    }


  }

  if (lambda[1] != 1)
  if (lambda[2] == 0 || lambda[3] == 0)
    warning("This function is not suitable for finding the c-optimal designs for estimating 'MED' or 'ED50' and the results may not be stable.")

  #####################################################################################################



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
  if (is.null(control$stop_rule))
    control$stop_rule <- "maxiter" # "equivalence" or "one_empire"
  if (is.null(control$stoptol))
    control$stoptol <- .99
  if (is.null(control$equivalence_every))
    control$equivalence_every <- 200
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



  ## inner space
  if (is.null(control$inner_space))
    control$inner_space <- "continuous" # or "vertices"
  if (!control$inner_space %in% c("continuous", "vertices", "discrete"))
    stop("inner space can be \"continuous\", \"vertices\" or \"discrete\"")
  if (control$inner_space == "discrete" && is.null(control$param_set))
    stop("'param_set' must be given")
  if (is.null(control$inner_maxeval))
    control$inner_maxeval <- 600
  if (is.null(control$check_inner_maxeval))
    control$check_inner_maxeval <- TRUE
  if (is.null(control$plot_cost))
    control$plot_cost <- FALSE
  if (is.null(control$plot_deriv))
    control$plot_deriv <- TRUE
  ##trace
  if (is.null(control$trace))
    control$trace <- FALSE


  ## not from the user
  if (is.null(control$maxeval_equivalence))
    control$maxeval_equivalence <- 6000 ## maxiter for directL that is used to find the maximum of the derivative function

  if (is.null(control$rseed))
    control$rseed <- NULL

  # not by user
  control$answering_merg_tol <- .005


  ### do not change the seed
  if (exists(".Random.seed")) {
    OldSeed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", OldSeed, envir = .GlobalEnv))
  }


  ##################################### warnings: fimfunc is not a character anymore
  ##  maybe you should change the fimfunc to fimchar later

  # fimfunc can be character or function. however for now it is only character.
  #if character, set fimchar and find the appropriate FIM function.
  if (is.character(fimfunc)) {
    fimchar <- fimfunc
    # change 'fimfunc' to be a function corresponding to 'fimmchar'.
    fimfunc <- check_lp_fimchar(fimchar = fimchar, lp = lp)
  }else
    fimchar <- "unknown"

  # to handle ...
  fimfunc2 <- function(x, w, param) {
    #we have no .. currently because the model is fixed!!
    out <- fimfunc(x = x, w = w, param = param)
    # out <- fimfunc(x = x, w = w, param = param,...)
    return(out)
  }

  # global variables needed for the definition of crfunc
  n_independent <- length(lx)
  npar <- length(lp)


  temp <- create_multiple (model = fimchar, fimfunc2 = fimfunc2,  type = type, multiple = multiple)
  crfunc <- temp$crfunc

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

  if (!is.null(initial)) {
    # convert the intial to matrix if it is a vector!
    if (!is.matrix(initial))
      initial <- t(as.matrix(initial))
    initial <- round(initial, 12) #because it may produce strange error when cheking the lower upper bound!!!
    # check if the length is true!
    if (!is.na(initial) && dim(initial)[2] != length(ld))
      stop("The number of columns of 'initial' does not match with length of countries and should be", length(ld))
    #  we use round to protect the function from strange behaviour
    initial_out <-  sapply(1:dim(initial)[1], function(j) any(round(initial[j,], 5) > round(ud, 5)) || any(round(initial[j,], 5) < round(ld, 5)))
    if (any(initial_out))
      stop("The initial vaule(s) in row ", paste(which(initial_out), collapse = ", "), " are (is) out of bound.")

  }

  ##### for multiobjective optimal design
  vrunif <- Vectorize(runif)
  endpoint <- matrix(runif(length(ld), ld, ud), 1, length(ld))
  endpoint[1, 1:2] <- c(lx, ux)
 if (is.null(initial))
   initial <- endpoint else
     initial <- rbind(initial, endpoint)

  ## the variables that will be added to control not by user, but by mica
  arg <-
    list(
      lx = lx, ux = ux, lp = lp, up = up, k = k, ld = ld, ud = ud,
      FIM = fimfunc2, crfunc = crfunc, locally = NULL,
      initial = initial,control = control, type = type,
      param = param, delta = delta, lambda = lambda
    )
  ## updating will be added to arg in iterate


  if (type == "multiple_locally"){
    arg$Psi_x <- temp$PsiMulti_x
    arg$Psi_Mu <- temp$PsiMulti_Mu
  }


  ICA_object <- list(arg = arg,
                     evol = NULL)
  class(ICA_object) <- c("list", "ICA")

  out <- iterate.ICA(object = ICA_object, iter = iter)
  return(out)

}

