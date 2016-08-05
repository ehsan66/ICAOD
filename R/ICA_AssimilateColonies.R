


##if  ExceedStrategy = Random, the points that exceed lower and upper bound take the random values between the bounds.
##if ExceedStrategy = Lower_Upper, the points that exceed lower nad upper bound take the lower and upper bounds.

AssimilateColonies2 <- function(TheEmpire,
                                AssimilationCoefficient,
                                VarMin,
                                VarMax,
                                ExceedStrategy,
                                AsssimilationStrategy,
                                sym,
                                MoveOnlyWhenImprove,
                                fixed_arg,
                                w_id,
                                equal_weight){
  #   TheEmpire: a list  contains the position of the colonies and imperialists and the cost values
  #   AssimilationCoefficient
  #   VarMin: lower bound of the positions
  #   VarMax: upper bounds of the positions
  #   ExceedStrategy: passed to 'CheckBoxConstraints' function. can be 'random', 'perturbed'
  #   AsssimilationStrategy: 'PICA' or original 'ICA' or 'OICA'
  #   sym: if the positions must be symetric
  #   MoveOnlyWhenImprove.
  #   fixed_arg. arguments that are fixed and will be passed to the optimization over the inner problem, passed to Calculate cost
  #   w_id: column index of weights in the position matrix. you can take w_id from fixed_arg, but is avoided to not make any confusion

  # Rerurn the assimilated positions of colonies and update  TheEmpire$ColoniesPosition, TheEmpire$ColoniesCost and   TheEmpire$ColoniesInnerParam


  NumOfColonies <- dim(TheEmpire$ColoniesPosition)[1]
  ##MARGIN is 2 since it should be subtracted columnwise!
  diff_vec <- -sweep(x = TheEmpire$ColoniesPosition, MARGIN = 2, STATS = TheEmpire$ImperialistPosition, FUN = "-")


  if(AsssimilationStrategy =="ICA")
    NewPosition = TheEmpire$ColoniesPosition +
    AssimilationCoefficient * matrix(runif(length(diff_vec)), dim(diff_vec)[1], dim(diff_vec)[2]) * diff_vec

  if(AsssimilationStrategy == "PICA")
    NewPosition = TheEmpire$ColoniesPosition +
    (AssimilationCoefficient * matrix(runif(length(diff_vec)), dim(diff_vec)[1], dim(diff_vec)[2]) - 1) * diff_vec


  NewPosition <- CheckBoxConstraints(Position = NewPosition, Lower = VarMin, Upper = VarMax, Strategy = "perturbed")


  ### for weights position
  if (!equal_weight){
    ##if you want to only move in feasible region!
    w_mat <- NewPosition[, w_id, drop = FALSE]
    even_odd <- ifelse(length(VarMin)%% 2 == 0, "even", "odd")
    NewPosition[, w_id] <- SumToOne(w_mat = w_mat, sym = sym, even_odd = even_odd)
  }

  temp <- Calculate_Cost(mat=NewPosition, fixed_arg=fixed_arg)
  nfeval <- temp$nfeval
  NewCost <- temp$cost
  nimprove <- 0 ## if there is no improvement this will be returnEd!!

  if(MoveOnlyWhenImprove){
    OldCost <- TheEmpire$ColoniesCost
    replace_id <- which(NewCost < OldCost)
    if(length(replace_id) != 0){
      TheEmpire$ColoniesPosition[replace_id,] <- NewPosition[replace_id, , drop = FALSE]
      TheEmpire$ColoniesCost[replace_id] <- NewCost[replace_id]
      TheEmpire$ColoniesInnerParam[replace_id, ] <- temp$inner_optima[replace_id, , drop = FALSE]
      nimprove <- length(replace_id)
    }##else dont change anything
  }
  if(!MoveOnlyWhenImprove){
    TheEmpire$ColoniesPosition <- NewPosition
    TheEmpire$ColoniesCost <- NewCost
    TheEmpire$ColoniesInnerParam <- temp$inner_optima
  }


  return(list(TheEmpire = TheEmpire, nfeval = nfeval, nimprove = nimprove))

}


