RevolveColonies <- function(TheEmpire, RevolutionRate, NumOfCountries,
                            lower,
                            upper,
                            sym,
                            sym_point,
                            fixed_arg,
                            w_id,
                            equal_weight){

  # Revolve the colonies and update 'TheEmpire'
  # the position is changed if and only if improvement happens

  nrevol <- 0 ## number of revolutions
  nfeval <- 0
  # the number of design points is odd or even
  even_odd <- ifelse(length(lower)%% 2 == 0, "even", "odd")
  sigma <- 0.1 * (upper - lower)
  n_col <- dim(TheEmpire$ImperialistPosition)[2]

  ## we choose the columns by random
  DimSize <- sample(x = 1:n_col, size = 1)
  DimIndex <- sample(x = 1:n_col, size = DimSize )
  new_imp <- TheEmpire$ImperialistPosition



  new_imp[, DimIndex]  <- new_imp[, DimIndex, drop = FALSE] + sigma[DimIndex] * runif(DimSize, -1, 1)

  new_imp[, DimIndex] <- new_imp[, DimIndex]  * (( new_imp [, DimIndex] <= upper[DimIndex]) & ( new_imp[, DimIndex] >= lower[DimIndex])) +
    (new_imp[, DimIndex] > upper[DimIndex]) * (upper[DimIndex] - .25 * (upper[DimIndex] - lower[DimIndex]) * runif(1)) +
    (new_imp[, DimIndex] < lower[DimIndex]) * (lower[DimIndex] + .25 * (upper[DimIndex] - lower[DimIndex]) * runif(1))

  if(!equal_weight)
    new_imp[, w_id] <- SumToOne(w_mat =  new_imp[, w_id, drop = FALSE], sym=sym, even_odd=even_odd)


  temp <- Calculate_Cost (mat = new_imp, fixed_arg = fixed_arg)
  nfeval <- nfeval + temp$nfeval
  new_imp_cost <- temp$cost
  #new_imp_cost <- -temp$cost
  if(new_imp_cost < TheEmpire$ImperialistCost){
    nrevol <- nrevol + 1
    TheEmpire$ImperialistCost <- new_imp_cost
    TheEmpire$ImperialistPosition <- new_imp
    TheEmpire$ImperialistInnerParam <- temp$inner_optima ## new_imp_inner_param
  }

  ##now revolve colonies
  NumOfRevolvingColonies <- round(RevolutionRate * length(TheEmpire$ColoniesCost))

  if(NumOfRevolvingColonies != 0){


    ##The colonies that should be revolved.
    RandomIndex <- sample(x = 1:length(TheEmpire$ColoniesCost),
                          size = length(TheEmpire$ColoniesCost),
                          replace = FALSE)[1:NumOfRevolvingColonies]


    ######################################################
    ### a local search for revolved colonies
    DimSize <- sample(x = 1:n_col, size = 1)
    DimIndex <- sample(x = 1:n_col, size = DimSize )

    new_colonies <- TheEmpire$ColoniesPosition[RandomIndex, , drop = FALSE]
    new_colonies[, DimIndex] <-  new_colonies[, DimIndex, drop = FALSE] +
      matrix(sigma[DimIndex], nrow = length(RandomIndex), ncol = DimSize, byrow = TRUE) * matrix(runif(length(new_colonies[, DimIndex, drop = FALSE])), length(  RandomIndex) , DimSize )
    if(!equal_weight)
      new_colonies[,w_id] <- SumToOne(w_mat = new_colonies[, w_id, drop = FALSE], sym=sym, even_odd=even_odd)
    ###check the lower upper bound

    new_colonies[, DimIndex]<- CheckBoxConstraints(Position = new_colonies[, DimIndex, drop = FALSE],
                                                   Lower = lower[DimIndex] , Upper = upper[DimIndex], Strategy = "perturbed")

    temp1 <- Calculate_Cost (mat = new_colonies, fixed_arg = fixed_arg)
    nfeval <- nfeval + temp$nfeval
    new_colonies_cost <- temp1$cost

    change_id <-  which(new_colonies_cost < TheEmpire$ColoniesCost[RandomIndex])
    if(length(change_id) != 0){
      nrevol <- nrevol + length(change_id)
      TheEmpire$ColoniesPosition[RandomIndex[change_id], ] <- new_colonies[change_id, , drop = FALSE]
      TheEmpire$ColoniesCost[RandomIndex[change_id]] <- new_colonies_cost[change_id]
      TheEmpire$ColoniesInnerParam[RandomIndex[change_id], ] <- temp1$inner_optima[change_id, , drop = FALSE]


    }
    ##################
  }
  return(list(TheEmpire = TheEmpire, nfeval = nfeval, nrevol = nrevol))
}
