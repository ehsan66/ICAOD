#' Imperialist Competitive Algorithm to find locally, minimax and standardized maximin D-optimal designs for nonlinear models
#'
#'
#'  Let \eqn{\Xi} be the space of all  approximate designs with
#'  \eqn{k} support points  at \eqn{x_1, x_2, ...,  x_k}{x1, x2, ...,  xk} from  design space \eqn{\chi} with
#'  corresponding weights  \eqn{w_1, . . . ,w_k}{w1, . . . ,wk}. Let \eqn{M(\xi, \theta)} be the Fisher information
#'  matrix (FIM) of a \eqn{k-}point design \eqn{\xi} and \eqn{\theta} be the vector of unknown parameters.
#'  A  minimax D-optimal design \eqn{\xi^*}{\xi*} is defined by
#'  \deqn{\xi^* = \arg \min_{\xi  \in \Xi}\max_{\theta \in \Theta} -\log|M(\xi, \theta)|.}{
#'  \xi* = arg min over \Xi {max over \Theta -log|M(\xi, \theta)|}.}
#'  Throughout, the terms ``inner problem'' and ``outer problem'' are used for optimization over \eqn{\Theta} and \eqn{\Xi}, respectively.\cr
#'  A standardized maximin D-optimal designs \eqn{\xi^*}{\xi*} is defined by
#'  \deqn{ \xi^* = \arg \max_{\xi \in \Xi}\inf_{\theta \in \Theta} \left[\left(\frac{|M(\xi, \theta)|}{|M(\xi_{\theta}, \theta)|}\right)^\frac{1}{p}\right],}{
#'   \xi* =  arg max over \Xi {inf over \Theta {|M(\xi, \theta)| / |M(\xi_\theta, \theta)|}^p},}
#'   where \eqn{p} is the number of model paramters and
#'    \eqn{\xi_{\theta}}{\xi_\theta}  is the locally D-optimal design with respect to \eqn{\theta}.\cr
#' A locally D-optimal designs \eqn{\xi^*}{\xi*} is defined by
#' \deqn{\xi^* = \arg \min_{\xi  \in \Xi} -\log|M(\xi, \theta_0)|.}{\xi* = arg min -log|M(\xi, \theta0)|,}
#' where the minization is over \eqn{\Xi} and \eqn{\theta_0}{\theta0} is the initial values for the unknown parameters.\cr
#'
#'
#'
#' @param fimfunc FIM as a \code{function} or \code{character}.
#' As a \code{character}, it must be the name of the FIM function that are available,
#'  e.g. \code{as.character(substitute(FIM_logistic))}.
#'    As a \code{function} it must return the information matrix.  See "Details".
#' @param lx lower bound of the design space \eqn{\chi}.
#' @param ux upper bound of the design space \eqn{\chi}.
#' @param lp lower bound of the region of unceratinty \eqn{\Theta}. The order of the dimension is the same as the order of the parameters in the argument \code{param}of \code{fimfunc}.
#' @param up upper bound of the region of unceratinty \eqn{\Theta}. If \code{lp = up}, then \eqn{\Theta = \theta_0}{\Theta = \theta0} and the genearted design is locally D-optimal design. See "Examples".
#' @param iter maximum number of iterations.
#' @param k number of design (support) points. Must be larger than the number of model parameters \eqn{p} to avoid singularity of the FIM.
#' @param control a list of control parameters. See "Details".
#' @param control_gosolnp  tuning parameters of function \code{\link[Rsolnp]{gosolnp}} for models that their locally optimal design do not have an
#'  analytical solution and is find by \code{\link[Rsolnp]{gosolnp}}. Only required when \code{type} is set to \code{'standardized'}.
#'  See "Details" of \code{\link{equivalence}}.
#' @param type a character strings; \code{"minimax"} for minimax optimal design, \code{"standardized"} for standardized maximin D-optimal design and \code{"locally"} for locally D-optimal design.
#' When \code{"locally"}, then  \code{lp} must be set equal to \code{up}.
#' @param initial a matrix of user intial countries or a vector of a country that will be inserted  into the initial countries of ICA. See "Details" .
#' @param locally locally a function that returns the value of determinant of FIM for
#'         the locally D-optimal design, i.e.  \eqn{|M(\xi_{\bold{\theta}}, \bold{\theta})|}{|M(\xi_\theta, \theta)|}.
#'          Only required when \code{type} is set to \code{"standardized"}. See "Details" of \code{\link{equivalence}}.
#' @param ... further arguments to be passed to the FIM function corresponding to \code{fimfunc}.
#'  For power logisitc model when \code{fimfunc} is equal to \code{"FIM_power_logistic"},
#'   the value of  \code{s} must be given here.
#' @details
#'
#'  \code{fimfunc}  as a \code{function} must have three arguments:
#'   1) desig points \code{x},
#'   2) weights \code{w} and 3) model parameters \code{param}.
#'    The output should be of type \code{matrix}.
#'     Further parameters can be set, but should be passed by  \code{...} in \code{mica} (like parameter \eqn{s} in power logistic model).
#'     See "Examples".\cr
#' \code{fimfunc} as a character string  must be the name of the FIM functions defined in this package.
#' The implemented FIM functions are given as follows:
#' \tabular{lll}{
#' \code{fiumfunc = "FIM_logistic"}            \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_logistic}} \cr
#' \code{fiumfunc = "FIM_logistic_4par"}       \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_logistic_4par}} \cr
#' \code{fiumfunc = "FIM_power_logistic"}      \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_power_logistic}} \cr
#' \code{fiumfunc = "FIM_michaelis"}           \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_michaelis}} \cr
#' \code{fiumfunc = "FIM_emax_3par"}           \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_emax_3par}} \cr
#' \code{fiumfunc = "FIM_loglin"}              \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_loglin}} \cr
#' \code{fiumfunc = "FIM_exp_2par"}            \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_exp_2par}} \cr
#' \code{fiumfunc = "FIM_exp_3par"}            \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_exp_3par}} \cr
#' \code{fiumfunc = "FIM_comp_inhibition"}     \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_comp_inhibition}} \cr
#' \code{fiumfunc = "FIM_noncomp_inhibition"}  \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_noncomp_inhibition}} \cr
#' \code{fiumfunc = "FIM_uncomp_inhibition"}   \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_uncomp_inhibition}} \cr
#' \code{fiumfunc = "FIM_mixed_inhibition"}    \tab  equivalent to    \tab \code{fiumfunc = \link{FIM_mixed_inhibition}} \cr
#' }
#' Setting \code{fimfunc} as a character strings only results in using the internal defined \code{locally} for standardized maximin
#' optimal designs. See each information matrix for information about the parameters and the analytical locally optimal design (if available).
#'
#'
#'
#'
#' The \code{control} argument is a list that can supply any of the following components:
#' \describe{
#'   \item{\code{ncount}}{number of countries.  Defaults to \code{40}.}
#'   \item{\code{nimp}}{number of imperialists. Defaults to 10 percent of \code{ncount}.}
#'   \item{\code{assim_coeff}}{assimilation coefficient. Defaults to \code{4}.}
#'   \item{\code{revol_rate}}{revolution rate. Defaults to \code{0.3}}
#'   \item{\code{damp}}{damp ratio for revolution rate. Less than one. \code{revol_rate} is decreased by \code{damp} in every iteration. Defaults to \code{0.99}}
#'   \item{\code{zeta}}{a coefficient to find the 'total cost' of empires. Defaults to \code{0.1}.}
#'   \item{\code{uniting_threshold}}{if the distance between two imperialists is less than the product of the \code{uniting_threshold} and the largest distance in the search space, then ICA unites these two empires.
#'   defaults to \code{0.02}.}
#'   \item{\code{assim_strategy}}{a character strings denotes assimilation strategy; \code{PICA} for perturbed ICA and \code{ICA} for the original version. Defaults to \code{PICA}.}
#'   \item{\code{lsearch}}{logical. Perform a local search on imperialists in every iteration?  Defaults to \code{TRUE}.}
#'   \item{\code{l}}{positive integer. the number of local search for each imperialist. Defaults to \code{2}.}
#'   \item{\code{only_improve}}{logical. In assimilation step, only move the colonies if the new position is better than the current position. Defaults to \code{TRUE}.}
#'   \item{\code{stop_rule}}{a character string denotes stopping rule.
#'   When \code{"maxiter"}, then  ICA only stops when  reachs the maximum number of iterations.
#'   When \code{"one_empire"}, then ICA stops if either all empires collapsed and one empire remains or reachs the maximum number of iterations.
#'   When \code{"equivalence"}, then ICA stops if either the D-efficiency lower bound (\code{DLB}) of the current design is greater than \code{stoptol} or reachs the maximum number of iterations.}
#'   \item{\code{stoptol}}{numeric between \eqn{0} and \eqn{1}. The minimum \code{DLB} for the best current imperialist (best design) to stop the algorithm by equivalence theorem when \code{stop_rule = "equivalence"}. Defaults to \code{0.99}.}
#'   \item{\code{equivalence_every}}{a positive integer. Check and compute \code{DLB} in every \code{equivalence_every} iteration. Checking equivalence theorem in small intervals slows down the algorithm. Defaults to \code{200}.}
#'   \item{\code{equal_weight}}{logical; whether the points should have equal weights. Defaults to \code{FALSE}.}
#'   \item{\code{sym}}{logical. Whether the design is symmetric around a point. If \code{TRUE} then \code{sym_point} must be given. Defaults to \code{FALSE}}
#'   \item{\code{sym_point}}{a vector of the same length as \code{lx}. The point that the design is symmetric around.
#'    Must have the same length of \code{lx}. See "Examples".}
#'   \item{\code{inner_space}}{a character string denote the inner space. Can be  \code{"continuous"}, \code{"vertices"} or \code{"discrete"}. See below. Defaults to \code{"continuous"}.}
#'   \item{\code{param_set}}{a matrix denotes the fixed values for the parameters when \code{inner_space = "discrete"}. Each row of the matrix is one set (vector) of the values of the parameters,
#'    i.e.  \eqn{\theta_{01}}{\theta01} in the first row,  \eqn{\theta_{02}}{\theta02} in the second row and so on. The number of columns should be equal to \code{length(lp)}. See "Examples".}
#'   \item{\code{inner_maxeval}}{maximum number of function evaluations for the continuous inner problem. It comes from the tuning parameters of  \code{\link[nloptr]{directL}}.
#'    Its value should be large enough to not miss any global optima. It is only applicable for standardized maxinim and minimax
#'    optimal designs. Defaults to \code{600}.}
#   \item{\code{check_inner_maxeval}}{logical. Increasing the current value of \code{inner_maxeval} in \code{\link[nloptr]{directL}} to solve the inner problem  given the output (best) design will change the value of criterion?
#   If yes a warning message will be printed; probably in the inner problem the value of \code{inner_maxeval}  was not large enough to assure finding the global optimum of the inner problem.
#   In this case, the algorithm should be re-runned with a larger value of \code{inner_maxeval}  that is given in the warning message or is saved in \code{best$msg}.}
#'   \item{\code{plot_deriv}}{logical. Should derivative be plotted whenever equivalence theorem is checked? Defaults to \code{TRUE}.}
#'   \item{\code{plot_cost}}{logical. Should the evolution of ICA, i.e. mean cost of all imperialists and cost of the best imperialist, be plotted in every iteration? Defaults to \code{TRUE}; when \code{type = "locally"} is \code{FALSE}.}
#'   \item{\code{trace}}{logical. Should the best generated design (best imperialist) and corresponding algorithm parameters be printed in every iteration? Defauls to \code{TRUE}.}
#'   \item{\code{n.seg}}{a positive integer required when checking the equivalence theorem to construct the answering set.
#'   The number of initial starting points for local optimzer (\code{\link[stats]{optim}}) to find all minima of the criterion on parameters space is equal to \code{(n.seg + 1)^length(lp)}. See "Details". Defaults to \code{4}}
#' }
#'
#'
#' Each row of \code{intial} is one design, i.e. concatenation of
#'  \eqn{\bold{x} = (x_1,...x_k)}{x = (x1,...xk)} and \eqn{\bold{w} = (w_1,...,w_k)}{w = (w1,...,wk)}.
#' The number of columns of \code{initial} is equal to \eqn{k \times n + k}{k times  n + k},
#'  where \eqn{n} is the number of model explanatory variables. For multi-dimensional design space, \eqn{\bold{x}} must be
#'  given dimension by dimension. See description of argument \code{x} in \code{\link{equivalence}}.\cr
# Please note that even when \code{equal_weight = TRUE},
# then \code{initial} must be of dimension \eqn{k \times n}{k times  n} and the columns corresponing
#  the weights must be given.
#'
#'
#'
#' In \code{control}, if \code{inner_space = "continuous"}, then the inner problem is an optimization over
#'  \eqn{\Theta =} (\code{lp}, \code{up}).
#'  If \code{inner_space = "discrete"}, then the inner problem is a discrete optimization problem over
#'   a set of initial values for the parameters, i.e.
#'    \eqn{\Theta = \{\theta_{01}, \theta_{02},...\}}{\Theta = {\theta01, \theta02, ....}}.
#'  In this case, the set of initial parameters should be given through \code{param_set}.
#'  If \code{inner_space = "vertices"} then the set of intial parameters are the vertices of \eqn{\Theta}. This should be set when the user is certain that the D-criterion is convex with respect to the parameter space for every given design.
#'  Please note that regardless of what \code{inner_space} is, checking the equivalence theorem  is done on continuous parameter space \eqn{\Theta}. See "Examples" on how to use this option. \cr
#'
#'
#'
#' For large parameter space or complex models it is important to increase \code{ncount},
#'  \code{inner_maxeval} and \code{n.seg} (for checking the equivalence theorem).\cr
#'
#' Please note that the speed of \code{mica} highly depends on the \code{inner_maxeval} parameter
#'  in \code{control} when \code{inner_space = "continuous"}.
#' \code{equivalence_every} and \code{l} in local search are the other factors that impact the CPU time.\cr
#'
#' From Section "Value",  note that \code{all_optima}, \code{all_optima_cost}, \code{answering}, \code{answering_cost},
#'    \code{mu},  \code{max_deriv} and \code{DLB} are  \code{NA} when the
#'       equivalence theorem was not requested for that iteration by \code{equivalence_every} in control.
#'        For example, if \code{equivalence_every = 100}, then  the equivalence theorem is only check for
#'        best designs in iteration \eqn{100}, \eqn{200}, \eqn{300} and so on.\cr
#' From Section "Value", \code{inner_param}  is equal to \eqn{\arg \max_{\theta \in \Theta} -\log|M(\xi, \theta)|,}{arg  max -log|M(\xi, \theta)|,} for minimax and
#'        \eqn{  \arg \inf_{\theta \in \Theta} \left[\left(\frac{|M(\xi, \theta)|}{|M(\xi_{\theta}, \theta)|}\right)^\frac{1}{p}\right],}{ =  arg inf {|M(\xi, \theta)| / |M(\xi_\theta, \theta)|}^p,}
#'        for standardized maximin D-optimal designs.
#'
#'
#'
#'
#'
#'
#Changing the \code{inner_space} to \code{"vertices"} or \code{"discrete"} can or requesting locally D-optimal design can
#  considerably  decrease the CPU time since the continuous optimization in the inner problem will be removed. But, as minimax approach is a conservative approach,
#  we recommend the user to be conservative as well and avoid making compromises to discretize \eqn{\Theta} or turning off the local search, check of the maximum number of evaluations in the inner problem or
#   using a very small value for \code{inner_maxeval} (unless the user is certain that the optima over \eqn{\Theta} are attained on the vertices, i.e. the criterion is convex with respect to \eqn{\Theta} for any given design).
#
#'
#' @return
#'  an object of class "ICA" that is a list contains another three lists:
#' \describe{
#'   \item{\code{arg}}{a list contains the arguments. Required for other methods.}
#'   \item{\code{evol}}{a list as length as the number of iterations to save the best design (best imperialist) in each iteration.
#'    \code{evol[[iter]]} contains:
#'     \tabular{lll}{
#'       \code{iter}                   \tab      \tab iteration. \cr
#'       \code{x}                      \tab      \tab design point. \cr
#'       \code{w}                      \tab      \tab design weight. \cr
#'       \code{min_cost}               \tab      \tab cost of the best imperialist. \cr
#'       \code{mean_cost}              \tab      \tab mean of costs of all imperialists. \cr
#'       \code{all_optima}             \tab      \tab all optima of the inner problem. \code{NA} for locally optimal design. \cr
#'       \code{all_optima_cost}        \tab      \tab cost of all optima of the inner problem. \code{NA} for locally optimal design. \cr
#'       \code{answering}              \tab      \tab answering set. \code{NA} for locally optimal design. \cr
#'       \code{answering_cost}         \tab      \tab cost of each element of answering set. \code{NA} for locally optimal design. \cr
#'       \code{mu}                     \tab      \tab found probability measure on answering set. \code{NA} for locally optimal design. \cr
#'       \code{max_deriv}              \tab      \tab maximum of the sensitivity function.  \cr
#'       \code{DLB}                    \tab      \tab D-efficiency lower bound. \cr
#'       \code{inner_param}            \tab      \tab inner parameter. See "Details".\cr
#'     }
#'   }
#'
#'   \item{\code{empires}}{a list of empires.}
#'   \item{\code{alg}}{a list that contains the best solution with components:
#'     \tabular{lll}{
#'       \code{nfeval}           \tab      \tab number of function evaluations. (the call for checking equivalence theorem is not counted) \cr
#'       \code{nlocal}           \tab      \tab number of successful local search. \cr
#'       \code{nrevol}           \tab      \tab number of successful revolutions. \cr
#'       \code{nimprove}         \tab      \tab number of successful movements toward the imperialists in assimilation step. \cr
#'       \code{convergence}      \tab      \tab Why the algorithm has stopped?. \cr
#'     }
#'   }
#' }
#' The design in \code{best} is the same best imperialist generated in last iteration and  stored in \code{evol}.
#'
#' @references
#' Masoudi, E., Holling, H., & Wong, W. K.  (in press). Application of imperialist competitive algorithm to find minimax and standardized maximin optimal designs. Computational Statistics & Data Analysis.
#' @examples
#'#######################################################################
#'## some examples for exponential model
#' \dontrun{
#'# finding standardized maximin D-optimal design
#'res <- mica(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(1, 1), up = c(1, 5),
#'     iter = 100, k = 3, type = "standardized", control = list(seed = 215))
#'res <- iterate(res, 10)
#'plot(res)
#'# finding minimax D-optimal design
#'mica(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(1, 1), up = c(1, 5),
#'     iter = 100, k = 3, type = "minimax", control = list(seed = 215))
#'}
#'# finding locally D-optimal design. Please note that  'lp' and 'up' are equal
#'mica(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(2, 3), up = c(2, 3),
#'     iter = 40, k = 2, type = "locally", control = list(seed = 215))
#'# locally D-optimal design is x1 = lx, x2 = 1/lp[2]
#'
#'# requesting an equally-weighted design, i.e w_1 = w_2 = ... w_k
#'res_loc <-mica(fimfunc = "FIM_exp_2par", lx = 0, ux = 1, lp = c(2, 3), up = c(2, 3),
#'               iter = 40, k = 2, type = "locally",
#'                control = list(seed = 215, equal_weight = TRUE))
#'\dontrun{
#'res_loc <- iterate(res_loc, 10) ## update the result
#'plot(res_loc)
#'# using symetric option for the logisitic model
#'mica(fimfunc = "FIM_logistic", lx = -5, ux = 5, lp = c(0, 1), up = c(3.5 , 1.25),
#'     iter = 100, k = 5, control = list(rseed = 215, sym = TRUE,
#'      sym_point = (0 + 3.5)/2),type = "minimax")
#'#######################################################################
#'# 2PL model
#'mica(fimfunc = "FIM_logistic", lx = -5, ux = 5, lp = c(-1, 1), up = c(1 , 2),
#'     iter = 100, k = 3, control = list(rseed = 215), type = "minimax")
#'
#'# an example on how to supply 'fimfunc' with a function
#'logistic_fim <- function(x, w, param){
#'  a <- param[1]
#'  b <- param[2]
#'  constant <- 1/(1 + exp(-b * (x - a)))
#'  constant <- constant * (1 - constant)
#'  A <-  sum(w * b^2 * constant)
#'  B <- sum(-w * b * (x - a) * constant)
#'  C <- sum(w * ((x -a)^2) * constant)
#'  mat <- matrix(c(A, B, B, C), 2, 2)
#'  return(mat)
#'}
#'
#'mica(fimfunc = logistic_fim, lx = -5, ux = 5, lp = c(-1, 1), up = c(1 , 2),
#'     iter = 100, k = 3,  control = list(rseed = 215), type = "minimax")
#'      ## is the same when 'fimfunc = "FIM_logistic'
#'
#'mica(fimfunc = logistic_fim, lx = -5, ux = 5, lp = c(-1, 1), up = c(-1 , 1),
#'     iter = 100, k = 3, control = list(rseed = 215), type = "locally")
#'#######################################################################
#'
#'#######################################################################
#'## how to use inner_space option in control list
#'
#'#### Enzyme kinetic models. examples for3D plots
#'mica(fimfunc = "FIM_comp_inhibition", lx = c(0, 0), ux = c(30, 60),
#'     lp = c(7, 4, 2), up = c(7, 5, 3), k =3, type = "standardized",
#'     iter = 300, control = list(rseed = 215, inner_maxit = 300,
#'                                stop_rule = "equivalence",
#'                                countries = 100, nimperialists = 10))
#'
#'## setting the parameter space as only the points on the vertices
#'mica(fimfunc = "FIM_comp_inhibition", lx = c(0, 0), ux = c(30, 60),
#'     lp = c(7, 4, 2), up = c(7, 5, 3), k =3, type = "standardized",
#'     iter = 300, control = list(rseed = 215, inner_space = "vertices",
#'                                stop_rule = "equivalence",
#'                                countries = 100, nimperialists = 10))
#'
#'## every row is one of the vertices of Theta
#'param_set <- matrix(c(7, 4, 2, 7, 5, 2, 7, 4, 3, 7, 5, 3),
#'                    ncol = 3, nrow = 4, byrow = TRUE)
#'res <-mica(fimfunc = "FIM_comp_inhibition", lx = c(0, 0), ux = c(30, 60),
#'           lp = c(7, 4, 2), up = c(7, 5, 3), k =3, type = "standardized",
#'           iter = 300, control = list(rseed = 215,inner_space = "discrete",
#'                                      stop_rule = "equivalence", countries = 100,
#'                                      nimperialists = 10, param_set = param_set))

#'
#'
#'#######################################################################
#'
#'#######################################################################
#'## optimal designs for the 1Pl model
## it has been proved that the locally optimal design is a
# one point design with weight 1
#'mica(fimfunc = "FIM_logistic_1par", lx = 0, ux = 5,
#'      lp = 2, up = 2, k = 3, iter = 100, type = "locally")
#'
#'lx <- -.5
#'ux <- .5
#'ux - lx <= 2 * log(2 + sqrt(3))
#'mica(fimfunc = "FIM_logistic_1par", lx = lx, ux = ux,
#'     lp = -1, up = 1, k = 1, iter = 10,
#'     type = "standardized")
#'
#'
#'
#'lx <- -2
#'ux <- 2
#'ux - lx <= 2 * log(2 + sqrt(3))
#'mica(fimfunc = "FIM_logistic_1par", lx = lx, ux = ux,
#'     lp = -1, up = 1, k = 1, iter = 10,
#'    type = "standardized")
#'#######################################################################
#'
#'}
#'
#' @export
####### full table
# \code{"logistic"}            \tab      \tab \code{\link{FIM_logistic}} \cr
# \code{"power_logistic"}      \tab      \tab \code{\link{FIM_power_logistic}} \cr
# \code{"michaelis"}           \tab      \tab \code{\link{FIM_michaelis}} \cr
# \code{"michaelis2"}          \tab      \tab \code{\link{FIM_michaelis2}} \cr
# \code{"michaelis3"}          \tab      \tab \code{\link{FIM_michaelis3}} \cr
# \code{"emax_3par"}           \tab      \tab \code{\link{FIM_emax_3par}} \cr
# \code{"poi1"}                \tab      \tab \code{\link{FIM_Poissson1}} \cr
# \code{"poi2"}                \tab      \tab \code{\link{FIM_Poissson2}} \cr
# \code{"second_poi"}          \tab      \tab \code{\link{FIM_second_poi}} \cr
# \code{"nbin"}                \tab      \tab \code{\link{FIM_nbin}} \cr
# \code{"loglin"}              \tab      \tab \code{\link{FIM_loglin}} \cr
# \code{"richards"}            \tab      \tab \code{\link{FIM_richards}} \cr
# \code{"weibull"}             \tab      \tab \code{\link{FIM_weibull}} \cr
# \code{"exp_2par"}            \tab      \tab \code{\link{FIM_exp_2par}} \cr
# \code{"exp_3par"}            \tab      \tab \code{\link{FIM_exp_3par}} \cr
# \code{"comp_inhibition"}     \tab      \tab \code{\link{FIM_comp_inhibition}} \cr
# \code{"noncomp_inhibition"}  \tab      \tab \code{\link{FIM_noncomp_inhibition}} \cr
# \code{"uncomp_inhibition"}   \tab      \tab \code{\link{FIM_uncomp_inhibition}} \cr
# \code{"mixed_inhibition"}    \tab      \tab \code{\link{FIM_mixed_inhibition}} \cr
# \code{"compartmental"}       \tab      \tab \code{\link{FIM_compartmental}} \cr
# \code{"sig_emax"}            \tab      \tab \code{\link{FIM_sig_emax}} \cr
# \code{"sig_emax1"}           \tab      \tab \code{\link{FIM_sig_emax1}} \cr
# \code{"sig_emax2"}           \tab      \tab \code{\link{FIM_sig_emax2}} \cr
# \code{"sig_emax3"}           \tab      \tab \code{\link{FIM_sig_emax3}} \cr
# \code{"sig_emax4"}           \tab      \tab \code{\link{FIM_sig_emax4}} \cr
# \code{"sig_emax5"}           \tab      \tab \code{\link{FIM_sig_emax5}} \cr
# \code{"sig_emax6"}           \tab      \tab \code{\link{FIM_sig_emax6}} \cr
# \code{"two_exp_censor1"}     \tab      \tab \code{\link{FIM_two_exp_censor1}} \cr
# \code{"two_exp_censor2"}     \tab      \tab \code{\link{FIM_two_exp_censor2}} \cr
# \code{"three_exp_censor1"}   \tab      \tab \code{\link{FIM_three_exp_censor1}} \cr
# \code{"three_exp_censor2"}   \tab      \tab \code{\link{FIM_three_exp_censor2}} \cr



mica <- function(fimfunc,
                 lx,
                 ux,
                 lp,
                 up,
                 iter,
                 k,
                 control = list(),
                 control_gosolnp = list(),
                 type,
                 initial = NULL,
                 locally = NULL,
                 ...) {

  #############################################################################
  #### common with mfw
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
  if (missing(iter))
    stop("\"iter\" is missing")
  if (missing(k))
    stop("\"k\" is missing")
  if (!is.function(fimfunc) && !is.character(fimfunc))
    stop(" \"fimfunc\" can be either \"character\" or \"function\"")
  if (length(lx) != length(ux))
    stop("Length of \"lx\" is not equal to length of \"ux\"")
  if (length(lp) != length(up))
    stop("length of \"lp\" is not equal to length of \"up\"")
  if (!is.numeric(k) || (k %% 1) != 0 || k <= 0)
    stop("\"k\" must be a positive integer number")
  if (k < length(lp))
    stop("\"k\" must be larger than the number of parameters to avoid singularity")
  if (!type %in% c("minimax", "standardized", "locally"))
    stop("\"type\" must be \"minimax\", \"standardized\" or  \"locally\"")
  if (!is.numeric(iter) || (iter %% 1) != 0 || iter <= 0)
    stop("\"iter\" must be a positive integer number")
  #############################################################################

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

  #############################################################################
  ## inner space
  if (is.null(control$inner_space))
    control$inner_space <- "continuous" # or "vertices"
  if (!control$inner_space %in% c("continuous", "vertices", "discrete"))
    stop("inner space can be \"continuous\", \"vertices\" or \"discrete\"")
  if (control$inner_space == "discrete" && is.null(control$param_set))
    stop("'param_set' must be given")
  if (is.null(control$inner_maxeval))
    control$inner_maxeval <- 600
  # if (is.null(control$check_inner_maxeval))
  #   control$check_inner_maxeval <- TRUE
  #############################################################################

  #############################################################################
  ## plots and trace and seed
  if (is.null(control$plot_cost)){
    control$plot_cost <- TRUE
    if (type == "locally")
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


  ##############################################################################
  ## if only one point design was requested, then the weight can only be one
  if (k == 1)
    control$equal_weight <- TRUE
  ##############################################################################


  #############################################################################
  # setting the Fisher information matrix
  #if character, set fimchar and find the appropriate FIM function.
  if (is.character(fimfunc)) {
    fimchar <- fimfunc
    fimfunc <- check_lp_fimchar(fimchar = fimchar, lp = lp)
  }else
    fimchar <- "unknown"

  # to handle ...
  fimfunc2 <- function(x, w, param) {
    out <- fimfunc(x = x, w = w, param = param,...)
    return(out)
  }
  #############################################################################

  #############################################################################
  ## some setting for 1PL (Rasch model)
  if (fimchar == "FIM_logistic_1par")
    control$inner_space <- "vertices"


  #############################################################################
  # create the locally form fimchar
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

  # global variables needed for the definition of crfunc
  n_independent <- length(lx)
  npar <- length(lp)

  #############################################################################
  ## creating the criterion function
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
    #if (!user_locally) ## is it necessary? why auxilary will not be used by user anyway!!
      auxiliary_locally <- list(lx = lx, ux = ux, npar =  npar, fimfunc = fimfunc, control = control_gosolnp)

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
        numerator <-
          det2(fimfunc2(
            x = x, w = w, param = param
          ), logarithm = FALSE)
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
  #############################################################################


  #############################################################################
  ## setting the criterion function
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
  #############################################################################

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
      lx = lx, ux = ux, lp = lp, up = up, k = k, ld = ld, ud = ud,
      FIM = fimfunc2, crfunc = crfunc, locally = locally,
      initial = initial,control = control, type = type
    )
  ## updating will be added to arg in iterate

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

