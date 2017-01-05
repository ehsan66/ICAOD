check_inner <- function(obj, maxeval){
  ### obj is an object of class "ICA"
  ## only when type = "minimax" or "standardized"
  it <- length(obj$evol)
 obj$evol[[it]]$x
 out <- nloptr::direct(fn = obj$arg$crfunc, lower = obj$arg$lp, upper = obj$arg$up,
                q = c(obj$evol[[it]]$x, obj$evol[[it]]$w),
                n_independent = length(obj$arg$lx),
                control=list(xtol_rel = .Machine$double.eps,
                             maxeval = maxeval) )
  return(out)
}

#
# check_inner(res, 500)$value ## 0.8257132
# check_inner(res, 1000)$value ## 0.8257132
# check_inner(res, 5000)$value ## 0.8257132
# check_inner(res, 10000)$value ## 0.8257132
# ## there is no change in the criterion value
