
Calculate_Cost <- function(mat, fixed_arg){

  # mat: is the matrix of positions
  # fixed_arg: passed to calculate



  if(!is.matrix(mat))
    stop("'mat' must be matrix")

  n_cost <- dim(mat)[1]
  if (fixed_arg$type != "optim_on_average")
  inner_optima_store <- matrix(NA, ncol = length(fixed_arg$lp), nrow=  n_cost) else
    inner_optima_store <- matrix(NA, ncol = dim(fixed_arg$param)[2], nrow=  n_cost)

  store <- vector("list", n_cost) #temporarily
  cost <- vector("double", n_cost)
  if(fixed_arg$equal_weight)
      w_equal <- rep(1/fixed_arg$k, fixed_arg$k)


  for(i in 1:dim(mat)[1]){
    x <- mat[i, fixed_arg$x_id]
    if(!fixed_arg$equal_weight)
    w <- mat[i, fixed_arg$w_id] else
      w <- w_equal
    if(fixed_arg$sym){

      x_w <- ICA_extract_x_w(x = x, w = w,
                             sym_point = fixed_arg$sym_point)
      x <- x_w$x
      w <- x_w$w
    }

    optim_func <- fixed_arg$optim_func
    store[[i]] <- optim_func(fn = fixed_arg$crfunc2,
                             x = x,
                             w = w,
                             lower = fixed_arg$lp,
                             upper = fixed_arg$up,
                             fixedpar = fixed_arg$fixedpar,
                             fixedpar_id = fixed_arg$fixedpar_id,
                             n_independent= fixed_arg$n_independent)
    MinRowId <- which.min(store[[i]]$minima[, fixed_arg$CostColumnId])
    cost[i] <- store[[i]]$minima[MinRowId , fixed_arg$CostColumnId]
    inner_optima_store[i, ] <- store[[i]]$minima[MinRowId , -fixed_arg$CostColumnId]

    x <- w <-MinRowId <- NA
  }

##we multiply the cost function by negative for maximin and minimax optimal design!

if(fixed_arg$type != "locally" && fixed_arg$type != "optim_on_average"){
  cost <- -cost
  #inner_optima_store <- matrix(NA, nrow = length(cost))
}


  nfeval <-  sum(sapply(store, "[[", "counts"))


  return(list(cost = cost, nfeval = nfeval, inner_optima = inner_optima_store))
}

