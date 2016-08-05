#' Checking minimax, standardized maximin and locally D-optimality of a givne design by equivalence theorem
#'
#'
#' An approximate design \eqn{\xi} is a probability measure defined on a user-selected design space \eqn{\chi}.
#' Let \eqn{\Xi} be the space of all such designs on \eqn{\chi} and let \eqn{\xi}
#'  be an approximate design with \eqn{k} support points  at \eqn{\boldsymbol{x}_1, \boldsymbol{x}_2, ..., \boldsymbol{x}_k}{x1, x2, ..., xk}
#'   from \eqn{\chi} with corresponding weights  \eqn{w_1, . . . ,w_k}{w1, ..., wk},
#'   \eqn{\sum_{i=1}^{k} w_i = 1}{w1 + w2 + ...+ wk = 1}.
#'    A design \eqn{\xi^*}{\xi*} is minimax D-optimal among all designs on \eqn{\chi} if and only if there exists a probability measure \eqn{\mu^*}{\mu*} on
#'    \deqn{A(\xi^*) = \left\{\nu \in \Theta \mid -log|M(\xi^*, \nu)| = \max_{\theta \in \Theta} -log|M(\xi^*, \theta)| \right\},}{
#'      A(\xi*) = {\nu belongs to \Theta | -log|M(\xi*, \nu)| = maxization with respect to \theta over \Theta on function -log|M(\xi*, \theta)|} ,}
#'        such that the following inequality holds for all \eqn{\boldsymbol{x} \in \chi}{x belong to \chi}
#'         \deqn{c(\boldsymbol{x}, \mu^*, \xi^*) = \int_{A(\xi^*)} tr M^{-1}(\xi^*, \nu)I(\boldsymbol{x}, \nu)\mu^* d(\nu)-p \leq 0,}{
#'          c(x, \mu*, \xi*) = integration over A(\xi*) with integrand  tr M^-1(\xi*, \nu)I(x, \nu)\mu* d(\nu)-p <= 0,}
#'           with equality at all support points of \eqn{\xi^*}{\xi*}.
#'            Here, \eqn{p} is the number of model parameters. \eqn{c(\boldsymbol{x}, \mu^*, \xi^*)}{c(x, \mu*, \xi*)} is called \bold{sensitivity function}.
#' The set \eqn{A(\xi^*)}{A(\xi*)} is sometimes called the \bold{answering set} of
#'  \eqn{\xi^*}{\xi*} and the measure \eqn{\mu^*}{\mu*} is a subgradient of the
#'    non-differentiable criterion evaluated at \eqn{M(\xi^*,\nu)}{M(\xi*,\nu)}.\cr
#'
#' @param fimfunc Fisher information matrix. Can be the name of the Fisher information matrix from FIM family functions available in this package as a
#'  character string or a function that returns the information matrix. See "Details" of \code{\link{mica}}.
#' @param x a vector of design points. When design space is multi-dimensional then \code{x} should be filled dimension by dimension. See "Examples".
#' @param w a vector of design weights.
#' @param lp lower bound of the region of unceratinty \eqn{\Theta}. The order of the dimension is the same as the order of the parameters in the argument \code{param}of \code{fimfunc}.
#' @param up upper bound of the region of unceratinty \eqn{\Theta}. If \code{lp = up}, then \eqn{\Theta = \boldsymbol{\theta_0} = lp}{\Theta = \theta_0 = lp} and the design is checked for locally D-optimality, See also argument \code{type}.
#' @param lx lower bound of the design space \eqn{\chi}.
#' @param ux upper bound of the design space \eqn{\chi}.
#' @param type a character strings; \code{'minimax'} for minimax optimal design, \code{'standardized'} for standardized maximin D-optimal design and \code{'locally'} for locally D-optimal design.
#' @param n.seg the number of starting points for finding all local optima to construct the answering set is \code{(n.seg + 1)^p}. Only for minimax and standardized maximin optimal designs. See "Details".
#' @param maxeval_equivalence maximum number of evaulations (\code{maxeval})  that will be passed to optimization function \code{\link[nloptr]{directL}} to find the maximum of the sensitivity function required for calculating DLB. See "Details".
#' @param locally a function that returns the value of determinant of FIM for
#'         the locally D-optimal design, i.e.  \eqn{|M(\xi_{\boldsymbol{\theta}}, \boldsymbol{\theta})|}{|M(\xi_\theta, \theta)|}.
#'          Only required when \code{type} is set to \code{"standardized"}. See "Details"
#' @param control_gosolnp  tuning parameters of function \code{\link[Rsolnp]{gosolnp}} for models that their locally optimal design do not have an
#'  analytical solution and is find by \code{\link[Rsolnp]{gosolnp}}. Only required when \code{type} is set to \code{'standardized'}. See "Details".
#' @param plot_3d a character strings to show which packages should be used for plotting the sensitivity function when
#' design space is of two dimension; \code{"lattice"} to use package \link[lattice]{lattice} and \code{"rgl"} to use package \link[rgl]{rgl}.
#' The package should be installed before use.
#' @param plot_sensitivity logical; sensitivity should be plotted? see "Details".
#' @param ... further argument to be passed to \code{fimfunc}.
#' @references
#' Atwood, C. L. (1969). Optimal and efficient designs of experiments. The Annals of Mathematical Statistics, 1570-1602.
#' @details
#'
#' For locally optimal designs the answering set has only one element that is \eqn{\nu = \boldsymbol{\theta_0} }{\nu = \theta_0} and \eqn{\mu =1}.
#' Thus, \code{n.seg} has no further use for locally optimal designs.\cr
#'
#'
#'     There is no theoretical rule on how to choose the number of points in \eqn{A(\xi^*)}{A(\xi*)} as support
#'  for the measure \eqn{\mu^*}{\mu*}  and they would have to be found by trial and error.
#'   To this end, we first find all the local maxima for optimization over \eqn{\Theta} by a local search (L-BFGS-B) with
#'    different \code{(n.seg + 1)^p} starting points  and then pick the ones nearest to the global minimum
#'     subject to a tolerance of \eqn{0.005}.\cr
#'
#'  If \eqn{\chi} is one or two dimensional, one may plot sensitivity function \eqn{c(\boldsymbol{x}, \xi^*, \mu^*)}{c(x, \xi*, \mu*)} versus
#'   \eqn{\boldsymbol{x} \in \chi}{x (belongs to \chi)} and  visually inspect whether the graph meets the conditions in the equivalence theorem.
#'     If it does, the design \eqn{\xi^*}{\xi*} is minimax optimal; otherwise it is not optimal.
#'
#'
#' We measure the closeness of a design
#'  \eqn{\xi} to the minimax optimal design  using its minimax D-efficiency defined by
#'   \deqn{d_{\mbox{eff}} = \left(\frac{\max_{\boldsymbol{\theta} \in \Theta} -\log|M(\xi_D, \boldsymbol{\theta})|}{\max_{\boldsymbol{\theta} \in \Theta} -\log|M(\xi, \boldsymbol{\theta})|}\right)^{1/p},}{
#'     d_eff = {(maximum over \Theta -log|M(\xi_D, \theta)|)/(maximum over \Theta -log|M(\xi, \theta)|)}^(1/p),}
#'       where \eqn{\xi_D} is  the minimax D-optimal design.
#' Using argument similar to Atwood (1969), we obtain a D-efficiency Lower Bound (DLB) for the minimax D-efficiency of a design \eqn{\xi} without knowing \eqn{\xi^*}{\xi*}.
#' DLB is caluclated by \eqn{p/(p + max_{\boldsymbol{x} \in \chi}c(\boldsymbol{x}, \mu, \xi))}{p/(p + maximum over \chi c(x, \mu, \xi))}, where \eqn{\mu^*}{\mu*}
#'     is the probability measure defined on \eqn{A(\xi)} that maximizes \eqn{c(x_\xi,\mu,\xi)}{c(x_\xi,\mu,\xi)} over all probability measures \eqn{\mu}
#'      and \eqn{x_\xi}{x_\xi} is any arbitrary support point of \eqn{\xi}.\cr
#'
#'  For standardazed maximin D-optimal design the answering set \eqn{N(\xi^*)}{N(\xi*)} is
#'    \deqn{N(\xi^*) = \left\{\boldsymbol{\nu} \in \Theta \mid \mbox{eff}_D(\xi^*, \boldsymbol{\nu}) = \min_{\boldsymbol{\theta} \in \Theta} \mbox{eff}_D(\xi^*, \boldsymbol{\theta}) \right\}.
#'      }{N(\xi*) = \{\nu belongs to \Theta  |eff_D(\xi*, \nu) = minimum over \Theta eff_D(\xi*, \theta) \},} where
#'      \eqn{\mbox{eff}_D(\xi, \boldsymbol{\theta}) =  (\frac{|M(\xi, \boldsymbol{\theta})|}{|M(\xi_{\boldsymbol{\theta}}, \boldsymbol{\theta})|})^\frac{1}{p}}{
#'      eff_D(\xi, \theta) =  (|M(\xi, \theta)|/|M(\xi_\theta, \theta)|)^(1/p)} and \eqn{\xi_\theta} is the locally D-optimal design with respect to \eqn{\theta}. \cr
#'  We measure the closeness of a design  \eqn{\xi} to the standardized maximin optimal design \eqn{\xi_D} using its standardized maximin  D-efficiency defined by
#'  \deqn{ d_{\mbox{eff}} = \frac{\min_{\boldsymbol{\theta} \in \Theta} \mbox{eff}_D(\xi, \boldsymbol{\theta})}{\min_{\boldsymbol{\theta} \in \Theta} \mbox{eff}_D(\xi_D, \boldsymbol{\theta})}.}{
#'  d_eff = (minimum over \Theta eff_D(\xi, \theta))/ (minimum over \Theta eff_D(\xi_D, \theta)).}
#'  Similar to the minimax design, we can also find standardized maximin D-efficiency
#'  lower bound for the generated standardized maximin design and locally D-optimal design.\cr
#'
#'  For standardized maximin D-optimal designs the function \code{locally} is created automatically when \code{locally = NULL}:
#'   \itemize{
#'      \item{when \code{fimfunc} is a \bold{character} strings from the available FIM and
#'       locally D-optimal design has a \bold{closed-form} for the chosen model: }{
#'        \code{locally} is an algebraic function that returns the value of determinant of the FIM for
#'         the locally D-optimal design, i.e.
#'         \eqn{|M(\xi_{\boldsymbol{\theta}}, \boldsymbol{\theta})|}{|M(\xi_\theta, \theta)|}. See "Details"
#'         of each defined FIM for the formula.}
#'      \item{when \code{fimfunc} is a \bold{character} strings but locally D-optimal design has \bold{no closed-form}:}{
#'         \code{\link[Rsolnp]{gosolnp}} is used to find the locally D-optimal design
#'          (within the class of minimally-supported and equally-weighted designs). \code{control_gosolnp} is a list of
#'          some of the tunning parameters of \code{\link[Rsolnp]{gosolnp}} that are: \code{n.sim}
#'           (default \code{500}), \code{n.restarts} (default \code{1})
#'            and \code{trace} (default \code{FALSE}).}
#'      \item{when \code{fimfunc} is a \bold{user given} function:}{ same as the previous case.}
#' }
#' User can also provide its own \code{locally}. In this case, \code{args(locally)}
#' must be \code{function (param, auxiliary)}, where \code{param} is the vector of parameters and
#' \code{auxiliary} is an obligatory arguments for internal use when \code{locally = NULL}. Please see
#' "Examples" on how to use \code{locally}.\cr
#'
#' For output \code{max_deriv}, if the local maximum is found (you can detect it from sensitivity plot) or DLB is negative,
#'  the value of \code{maxeval_equivalence} should be increased to return the global maximum.\cr
#'
#'  The criterion value for locally D-optimal design is
#'  \deqn{-\log|M(\xi, \boldsymbol{\theta_0} )|;}{-log|M(\xi, \theta_0 )|;} for minimax optimal design is
#'  \deqn{\max_{\theta \in \Theta} -\log|M(\xi, \theta)|;}{max -log|M(\xi, \theta)|;}
#'  for standardized maximin optimal design is
#'  \deqn{\inf_{\theta \in \Theta} \left[\left(\frac{|M(\xi, \theta)|}{|M(\xi_{\theta}, \theta)|}\right)^\frac{1}{p}\right].}{
#'   inf {|M(\xi, \theta)| / |M(\xi_\theta, \theta)|}^p.}
#'
#'
#'
#'
#' @examples
#'#############################################################
#'## check locally optimality: lp = up and type = "locally"
#'inipar <- c(2, 3)
#'equivalence (fimfunc = "FIM_logistic",
#'             x = c(1.485526, 2.51446 ),
#'             w = c(.5, .5),
#'             lx = -5, ux = 5,
#'             lp = inipar, up = inipar,
#'             type = "locally")
#'
#'##############################################################################
#'## standardized maximin D-optimal design does not depend on theta0 and theta1,
#'##  so we fix them locally D-optimal design has a closed-form which is defined
#'##  internally
#'equivalence (fimfunc = "FIM_loglin",
#'             x = c(0, 4.2494, 17.0324, 149.9090),
#'             w = c(0.3204, 0.1207, 0.2293, 0.3296),
#'             lx = 0, ux = 150,
#'             lp = c(2, 2, 1), up = c(2, 2, 15),
#'             type = "standardized")
#'
#'##############################################################################
#'\dontrun{
#'## there is no analytical solution for locally optimal design for this model
#'## gosolnp automatically will be used to find the locally
#'##  optimal design in the denominator of standardized criterion.
#'##  Becasue, it is two-level nested optimization
#'##  (first level on parameter space) and second level on design space)
#'##  it takes so long to find 'all_optima' and construct 'answerign' set.
#'equivalence (fimfunc = "FIM_power_logistic",
#'             x = c(-4.5515, 0.2130, 2.8075),
#'             w = c(0.4100, 0.3723, 0.2177),
#'             lx = -5, ux = 5,
#'             lp = c(0, 1), up = c(3, 1.5),
#'             type = "standardized",
#'             s = .2)
#' }
#'
#'
#'############################################################################
#'### when a design point is of two dimension
#'\dontrun{
#'equivalence (fimfunc = "FIM_mixed_inhibition",
#'             x = c(3.4614, 4.2801, 30, 30, 0, 3.1426, 0, 4.0373 ),
#'             w = rep(1/4, 4),
#'             lx = c(0, 0), ux = c(30, 60),
#'             lp = c(7, 4, 2, 4), up = c(7, 5, 3, 5),
#'             type = "standardized")
#'## here the design points are x1 = c(3.4614, 0), x2 = c(4.2801, 3.1426),
#'## x3 = c(30, 0), x4 = c(30, 4.0373)
#'## using package rgl (rgl must be installed for plot)
#'equivalence (fimfunc = "FIM_mixed_inhibition",
#'             x = c(3.4614, 4.2801, 30, 30, 0, 3.1426, 0, 4.0373 ),
#'             w = rep(1/4, 4),
#'             lx = c(0, 0), ux = c(30, 60),
#'             lp = c(7, 4, 2, 4), up = c(7, 5, 3, 5),
#'             type = "standardized", plot_3d = "rgl")
#'
#'equivalence (fimfunc = "FIM_comp_inhibition",
#'             x = c(3.4432, 30, 30, 0, 0, 18.8954),
#'             w = rep(1/3, 3),
#'             lx = c(0, 0), ux = c(30, 60),
#'             lp = c(7, 4, 2), up = c(7, 5, 3),
#'             type = "standardized")
#'}
#'##########################################################################
#'##########################################################################
#'## defining function 'locally'
#'locally_det<- function(param,  auxiliary){
#'  ## param is the vector of theta = (theta0, theta1, theta2)
#'  ux <- 0
#'  lx <- 150
#'  xstar <- (ux + param[3]) * (lx + param[3]) * (log(ux + param[3]) -
#'   log(lx + param[3]))/(ux - lx) - param[3]
#'  denominator <- det(FIM_loglin(x = c(lx, xstar, ux) , w = rep(1/3, 3), param = param))
#'  return(denominator)
#'}
#'equivalence (fimfunc = "FIM_loglin",
#'             x = c(0, 4.2494, 17.0324, 149.9090),
#'             w = c(0.3204, 0.1207, 0.2293, 0.3296),
#'             lx = 0, ux = 150,
#'             lp = c(2, 2, 1), up = c(2, 2, 15),
#'             locally = locally_det,
#'             type = "standardized")
#' @return
#'  an object of class \code{'equivalence'} that is a list contains:
#'  \describe{
#'  \item{\code{type}}{argumet \code{type} required for print methods.}
#'  \item{\code{all_optima}}{a matrix; all optima of the inner problem (optimization over the parameter space for the given design to find the minimum efficiency for standardized maximin design or maximum inefficiency for minimax. \code{NA} for locally optimal design.}
#'  \item{\code{all_optima_cost}}{cost of each element of \code{all_optima}. \code{NA} for locally optimal design.}
#'  \item{\code{answering}}{a matrix; answering set chosen from \code{all_optima}. \code{NA} for locally optimal design.}
#'  \item{\code{answering_cost}}{cost of each element of answering set. \code{NA} for locally optimal design.}
#'  \item{\code{mu}}{probability measure on answering set. Equal to \eqn{1} for locally D-optimal design.}
#'  \item{\code{max_deriv}}{maximum of the sensitivity function}
#'  \item{\code{DLB}}{D-efficiency lower bound. If negative, the value of \code{maxeval_equivalence} should be increased to find the global maximum.}
#'  \item{\code{crtval}}{criterion value. See "Details".}
#'  }
#' @export
#' @seealso \code{\link{print.equivalence}}, \code{\link{equivalence_on_average}} and \code{\link{equivalence_multiple}}.




equivalence <- function(fimfunc,
                        x, w,
                        lx, ux,
                        lp, up,
                        type,
                        n.seg = 6,
                        maxeval_equivalence = 6000,
                        locally = NULL,
                        control_gosolnp = NULL,
                        plot_3d = "lattice",
                        plot_sensitivity = TRUE,
                        ...){

  if (!is.function(fimfunc) && !is.character(fimfunc))
    stop("'fimfunc' can be either character or function")



  if (!type %in% c("minimax", "standardized", "locally"))
    stop("\"type\" must be \"minimax\", \"standardized\",  \"locally\"")

  if (length(lx) != length(ux))
    stop("Length of \"lx\" is not equal to length of \"ux\"")
  if (length(lp) != length(up))
    stop("length of \"lp\" is not equal to length of \"up\"")


  if (type == "standardized"){
    if (is.null(control_gosolnp$trace))
      control_gosolnp$trace <- FALSE
    if (is.null(control_gosolnp$n.restarts))
      control_gosolnp$n.restarts <- 1
    if (is.null(control_gosolnp$n.sim))
      control_gosolnp$n.sim <- 500
    if (is.null(control_gosolnp$rseed))
      control_gosolnp$rseed <- NULL
  }



  #############################################################################
  # fimchar

  # fimfunc can be character or function. however for now it is only character.
  #if character, set fimchar and find the appropriate FIM function.
  if (is.character(fimfunc)){
    fimchar <- fimfunc
    # change 'fimfunc' to be a function corresponding to 'fimmchar'.
    fimfunc <- check_lp_fimchar(fimchar = fimchar, lp = lp)
  }else
    fimchar <- "all"
  #############################################################################

  #############################################################################
  if (type == "standardized"){
    if(is.null(locally)){ ## if locally does not exist! we create the locally main function
      temp1 <- create_locally(fimchar = fimchar, lx = lx, ux = ux)
      locally <- temp1$locally_det
      user_locally <- FALSE
    } else {
      user_locally <- TRUE ## is the locally given by the user?
    }
    if (user_locally)
      error_msg <- "Please check your given 'locally' function." else {
        if (temp1$use_gosolnp)
          error_msg <- "please increase 'n.sim' and 'n.restarts' in 'control_gosolnp' list." else
            error_msg <- "Please contact the maintainer."
      }
  }
  #############################################################################


  if (type == "standardized" && is.null(locally))
    stop("Please specify the 'locally'")


  n_independent <- length(lx)
  npar <- length(lp)
  answering_merg_tol <- .005
  if (type == "locally")
    param_locally <- up



  fimfunc2 <- function(x, w, param){
    fimfunc(x, w, param, ...)
  }



  #############################################################################
  # definign criterion, coppied from 'mica'
  if (type == "minimax") {
    crfunc_minimax <- function(param, q, n_independent) {
      lq <- length(q) # q is the design points and weights
      pieces <- lq / (n_independent + 1)
      x_ind <- 1:(n_independent * pieces)
      w_ind <- (x_ind[length(x_ind)] + 1):lq
      x <- q[x_ind]
      w <- q[w_ind]

      minimax_crfunc <-
        det2(fimfunc2(x = x, w = w, param = param), logarithm = TRUE) - 5000 *
        (sum(w) - 1) ^ 2   ## -(-det+pen) = det-pen
      return(minimax_crfunc)
      locally = NULL
    }
    # we dont need minus for minimax because we want to maximize the -log(det) or minimze log(det)
  }
  if (type == "standardized") {
    auxiliary_locally <- list(lx = lx, ux = ux, npar =  npar, fimfunc = fimfunc2, control = control_gosolnp)
    crfunc_standardized <- function(q, param, n_independent) {
      lq <- length(q)
      pieces <- lq / (n_independent + 1)
      x_ind <- 1:(n_independent * pieces)
      w_ind <- (x_ind[length(x_ind)] + 1):lq
      x <- q[x_ind]
      w <- q[w_ind]
      # here we start a while loop to garantee that the locally optimal design is given by denominator.
      # we repeat finding the locally three times with different seeds and then we produce an error if for all three times the given design was not optimnal!
      continue <- TRUE
      counter <- 1
      while (continue) {
        denominator <- locally(
          param = param,
          auxiliary = auxiliary_locally
        )
        numerator <- det2(fimfunc2(x = x, w = w, param = param), logarithm = FALSE)
        fraction <- (numerator / denominator)
        if (round(fraction, 7) > 1) {
          #continue is still TRUE so the loop will be continuing
          counter <- counter + 1
        }else
          continue <- FALSE

        if (counter == 5) {
          stop("D-efficiency of the non-optimal design is higher than the optimal design for ",
               paste("theta_0 = c(", paste(round(param, 5), collapse = ","), ").",sep = ""),
               "\nProbably the generated design in 'locally' is not true locally D-optimal design with respect to theta_0. ",
               error_msg)
        }

      }
      if (npar %% 2 != 0) {
        maximin_crfunc <- (fraction) ^ (1 / npar)
      }else{
        maximin_crfunc <-  ifelse(fraction < 0, 0,(fraction) ^ (1 / npar))
      }
      return(maximin_crfunc)
    }
  }
  if (type == "locally") {
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
  }
  # set the criterion function
  if (type == "minimax")
    crfunc <- crfunc_minimax
  else
    if (type == "standardized")
      crfunc <- crfunc_standardized
  else
    if (type == "locally")
      crfunc <- crfunc_locally
  else
    stop("wrong 'type' is selected. Please check 'type'")
  ## end fo creating crfunc
  ##############################################################################



  ##############################################################################
  # dealing with the case when lp[i] = up[i]
  # coppied from iterate
  if (type != "locally"){
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
  ##############################################################################



  ######################################################################################################
  # required for finding the answering set for verification
  # copied from iterate
  #if (length(lp) <= 2)
  optim_starting <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  n_independent){
    out <- optim2(fn = fn, lower = lower, upper = upper,
                  n.seg = n.seg,
                  q = c(x, w),
                  fixedpar = fixedpar, fixedpar_id = fixedpar_id,
                  n_independent= n_independent)
    minima <- out$minima
    counts <- out$counts
    return(list(minima =minima, counts = counts))
  }
  ##############################################################################




  ## to just be in line with mica and do not mix Psi_mu and PSi_x
  arg <-
    list(
      lx = lx, ux = ux, lp = lp, up = up,
      FIM = fimfunc2)

  ##############################################################################
  ## Psi_mu and Psi_x{
  arg$Psi_x <- Psi_x
  arg$Psi_Mu <- Psi_Mu
  if (length(lx) == 2)
    arg$Psi_xy <- Psi_xy
  ##  Psi_x works for both one, two and three dimensional
  ## but Psi_xy is needed for plotting becasue the function should have two arguments

  ################################################################################
  ### Psi as a function of x and x, y for plotting. Psi_x defined as minus psi to find the minimum
  ## Psi_x is mult-dimensional, x can be of two dimesnion.
  Psi_x_minus <- function(x1, mu, FIM,  x, w,  answering){
    Out <- arg$Psi_x(x1 = x1, mu = mu, FIM = FIM,  x = x, w = w, answering = answering)
    return(-Out)
  }
  if(length(arg$lx) == 1)
    Psi_x_plot <-  arg$Psi_x ## for PlotPsi_x

  # it is necessary to distniguish between Psi_x for plotiing and finding DLB becasue in plotting for models with two
  # explanatory variables the function should be defined as a function of x, y (x, y here are the ploints to be plotted)

  if(length(arg$lx) == 2)
    Psi_x_plot <- arg$Psi_xy
  #when length(lx) == 1, then Psi_x_plot = Psi_x
  ################################################################################


  ##############################################################################
  # finding the answering set, measure and DLB
  # coppied from iterate
  if (type != "locally"){
    ## finding the answering set, measure

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
                             tol = answering_merg_tol,
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
    answering <- matrix(param_locally, nrow = 1) ## we need it for find measure and check. they use answering set
    mu <- 1 # we need it for check and plot
    answering_cost <- all_optima <- all_optima_cost <- NA
  }

  ##########################################################################
  # find the maximum of derivative function
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
                                          maxeval = maxeval_equivalence))

  ##sometimes the optimization can not detect maximum  in the bound, so here we add the cost values on the bound
  vertices_outer <- make_vertices(lower = arg$lx, upper = arg$ux)
  check_vertices <- find_on_points(fn = arg$Psi_x,
                                   points = vertices_outer,
                                   mu = mu,
                                   answering = answering,
                                   x = x,
                                   w = w,
                                   FIM = arg$FIM)
  check_vertices <- check_vertices$minima[, dim(check_vertices$minima)[2]]
  ## minus because for optimality check we minimize the minus derivative function
  max_deriv <- c(-OptimalityCheck$value, check_vertices)
  max_deriv <- max(max_deriv)
  ##########################################################################

  # D-efficiency lower bound
  DLB <- npar/(npar + max_deriv)


  if (plot_sensitivity)
    PlotPsi_x(lower = lx, upper =   ux,
              Psi_x = Psi_x_plot,
              FIM  = fimfunc2,
              mu = mu,
              x = x,
              w = w,
              plot_3d =plot_3d,
              answering = answering)



  object <- list(type = type,
                 all_optima = all_optima,
                 all_optima_cost = all_optima_cost,
                 answering = answering,
                 answering_cost = answering_cost,
                 mu = mu,
                 max_deriv = max_deriv,
                 DLB = DLB)


  if (type == "locally")
    object$crtval <- crfunc(param = lp, q = c(x, w), n_independent = 1)
  if (type == "minimax")
    object$crtval <- max(answering_cost)
  if(type == "standardized")
    object$crtval <- min(answering_cost)

  class(object) <- c("list", "equivalence")
  return(object)
}
