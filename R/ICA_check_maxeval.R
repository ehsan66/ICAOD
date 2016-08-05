## this function check whether a higher value of maxeval will give a better value or not
check_maxeval <- function(fn, lower, upper, maxeval, ...){
  if (missing(fn))
    stop("'fn' is missing")
  if (missing(fn))
    stop("'lower' is missing")
  if (missing(upper))
    stop("'upper' is missing")
  if (missing(maxeval))
    stop("'maxeval' is missing")
  ## WARNINGS: fn should have argument 'param'
  fn1 <- function(param) fn(param, ...)

  # fn1(param = c(1, 1))
  opt <-list()

  maxeval_vec <- maxeval + c(100, 200, 300, 500, 600, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 10000, 20000, 30000)
  for(i in 1:length(maxeval_vec))
    opt[[i]] <- directL(fn = fn1, lower = lower, upper = upper,
                       nl.info = FALSE,
                       control=list(xtol_rel=.Machine$double.eps,
                                    maxeval = maxeval_vec[i]))
  val <- sapply(opt, "[[", 'value')
  

  

  ## all equal?
  all.eq <- (abs(max(val) - min(val)) < .Machine$double.eps ^ 0.5)
  
  ############################################
  ## we produce a warnings if the points that directL missed is NOT on the vertices, because if it was on any of the vertices 
  ## cause no problem as our optimization function cosider this case as well.
  vertices <- make_vertices(lower = lower, upper = upper)
  vertices_val <- find_on_points(fn = fn1, points = vertices)
  vertices_val <- min(vertices_val$minima[, dim(vertices_val$minima)[2]])
  
  
  
  is.min_on_vertices <- ((vertices_val - min(val)) <  .Machine$double.eps ^ 0.5)
  ##################################################
  
  if (!all.eq && !is.min_on_vertices){
    temp_id <- which(abs(val - val[length(val)]) < .Machine$double.eps ^ 0.5)
    maxeval_recom <- maxeval_vec[temp_id[1]]
    msg <- paste("'inner_maxeval' recommended to be larger than ",   maxeval_recom, "\nPlease also consider increasing 'n.seg' if the value of 'inner_maxeval' is large.")
    warning(msg, call. = FALSE)
  }else{
    maxeval_recom <- maxeval
    msg <- NULL
  }
  return(list(msg = msg, maxeval = maxeval_recom, opt = opt))
}




#
# x <- c( 0.005676761, 0.02878325, 0.08840083, 0.2416213, 0.5159651, 0.7141474, 0.7154587, 0.8143749, 0.8831955, 0.9434678 )
# w <- c(0.2161778, 0.1500531, 0.1684801, 0.2107632, 0.09221182, 0.005554795, 0.002499739, 0.00726202, 0.08737721, 0.05962014)
# n_independent <- 1
# lower <- c(1, 1)
# upper <- c(1, 200)
#
#
#
# locally <- function(param, lx, ux, npar, fimfunc){
#
#   # warning: in the original paper the locally optimal design has only been found for lx = 0, but here by numerical results we
#   # saw that it is true for all lx
#   denominator <- det2(fimfunc(x = c(lx, 1/param[2]) , w = c(.5, .5), param = param) ,logarithm = FALSE)
#
#   return(denominator)
#
# }
# fn <- function(q, param, n_independent) {
#   lq <- length(q)
#   pieces <- lq / (n_independent + 1)
#   x_ind <- 1:(n_independent * pieces)
#   w_ind <- (x_ind[length(x_ind)] + 1):lq
#   x <- q[x_ind]
#   w <- q[w_ind]
#   npar = 2
#   # here we start a while loop to garantee that the locally optimal design is given by denominator.
#   # we repeat finding the locally three times with different seeds and then we produce an error if for all three times the given design was not optimnal!
#   continue <- TRUE
#   counter <- 1
#   while (continue) {
#     denominator <- locally(
#       param = param,
#       lx =0,
#       ux = 1,
#       npar =  2,
#       fimfunc = FIM_exp_2par
#     )
#     numerator <-
#       det2(FIM_exp_2par(
#         x = x, w = w, param = param
#       ), logarithm = FALSE)
#     fraction <- (numerator / denominator)
#     if (round(fraction, 7) > 1) {
#       #continue is still TRUE so the loop will be continuing
#       counter <- counter + 1
#     }else
#       continue <- FALSE
#
#     if (counter == 5) {
#       stop("\nD-efficiency of the non-optimal design is higher than the optimal design")
#       # "\nProbably the generated design in 'locally' is not true locally D-optimal design.")
#       #"\nIf the 'locally' was set by you, please check it and make sure that the returned design is optimal with respect to initial values of parameters",
#       #paste(" c(", paste(round(param, 5), collapse = ","), ").",sep = ""),
#       #"\nOtherwise, increase 'n.sim' or 'n.restarts' in 'gosolnp' control list.")
#     }
#
#   }
#   if (npar %% 2 != 0) {
#     maximin_crfunc <- (fraction) ^ (1 / npar)
#   }else{
#     maximin_crfunc <-  ifelse(fraction < 0, 0,(fraction) ^ (1 / npar))
#   }
#   return(maximin_crfunc)
# }
#
#
# test <- check_maxeval(fn, lower, upper, maxeval = 50, q = c(x, w), n_independent = 1)
# test
