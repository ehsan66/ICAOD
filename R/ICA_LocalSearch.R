##the LocalSearch is only for Point
LocalSearch <- function(TheEmpire, lower, upper, l, fixed_arg = fixed_arg){
  # the local search is done for half of the imperialists for 'l' times

  nfeval <- 0
  n_success <- 0
  imperialist <- as.vector(TheEmpire$ImperialistPosition)

  for(k in 1:(length(imperialist)/2)){

    counter <- 1
    while(counter <= l){
      diff_vec <- -sweep(x = TheEmpire$ColoniesPosition, MARGIN = 2, STATS = TheEmpire$ImperialistPosition, FUN = "-")

      ## 'd' is the maximal range is set to the distance between the imperialist and its closest colony in
      # the same emipire divided by  the square root of the number of veriables
      d <- apply(diff_vec, 1, function(x)sqrt(sum(x^2)))
      d <- min(d)/sqrt(length(imperialist)/2) #because the local search is only done for the half of the dimension

      if(round(d, 8) == 0)
        d <- .05 ## we set the d be equal to .05 if d is eqaul to zero!
      ## because if d is equal to zero then local search is not useful anymore!!!

        NewPos <-  TheEmpire$ImperialistPosition
        lambda <- runif(1, -1, 1)
        NewPos[k] <-  NewPos[k] + lambda * d * .05

        NewPos[k] <- (NewPos[k] <= upper[k] & NewPos[k] >= lower[k]) * NewPos[k] +
          (NewPos[k] > upper[k]) * (upper[k] - .25 * (upper[k] - lower[k]) * runif(1)) +
          (NewPos[k] < lower[k]) * (lower[k] + .25 * (upper[k] - lower[k]) * runif(1))


      output <- Calculate_Cost(mat = matrix(NewPos, nrow = 1), fixed_arg = fixed_arg)

      NewPos_cost <- output$cost
      fneval_candidate <- output$nfeval

      if( NewPos_cost < TheEmpire$ImperialistCost){
        n_success <- n_success + 1
        nfeval <- nfeval + fneval_candidate
        TheEmpire$ImperialistPosition <- matrix(NewPos, 1, length(lower))
        TheEmpire$ImperialistCost <- NewPos_cost
        TheEmpire$ImperialistInnerParam <- output$inner_optima
        counter <- l
      }
      counter <- counter + 1
    }
  }
  return(list(TheEmpire = TheEmpire, nfeval = nfeval, n_success = n_success))

}
