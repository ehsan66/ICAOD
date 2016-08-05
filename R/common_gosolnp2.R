#' @importFrom Rsolnp gosolnp
## gosolnp change the random seed or add.Random.seed when you use it. So we modified it.
gosolnp2 <- function(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL,
                     ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB = NULL,
                     UB = NULL, control = list(), distr = rep(1, length(LB)),
                     distr.opt = list(), n.restarts = 1, n.sim = 20000, cluster = NULL,
                     rseed = NULL, ...){
  if(exists(".Random.seed")){
    OldSeed <- get(".Random.seed", envir = .GlobalEnv)
    #if you call directly from update and not minimax!
    on.exit(assign(".Random.seed", OldSeed, envir = .GlobalEnv))

  }else
    on.exit(rm(.Random.seed, envir = .GlobalEnv)) ## becuse it produces .Random.seed

# becuase gosolnp uses the system.time() to set the rseed if it is NULL
# and it causes the results in while loop in minimax_FW for finding locally optimal design in maximin type be the same because of
# the rseed.
  if(is.null(rseed))
    rseed <- runif(1)

  out <- gosolnp(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun, eqB = eqB,
                 ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB, LB = LB,
                 UB = UB, control = control, distr = distr,
                 distr.opt = distr.opt, n.restarts = n.restarts, n.sim = n.sim, cluster = cluster,
                 rseed = rseed,...)

  return(out)
}

