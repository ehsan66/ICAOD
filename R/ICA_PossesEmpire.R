PossesEmpire <- function(TheEmpire){

  ## Exchanging Position of the Imperialist and  the best Colony, if the best colony is better

  ColoniesCost <- TheEmpire$ColoniesCost

  BestColonyInd <- which.min(ColoniesCost)

  if(ColoniesCost[BestColonyInd] < TheEmpire$ImperialistCost){

  OldImperialistPosition <- TheEmpire$ImperialistPosition
  OldImperialistCost <- TheEmpire$ImperialistCost
  OldImperialistInnerParam <- TheEmpire$ImperialistInnerParam

  TheEmpire$ImperialistPosition <- TheEmpire$ColoniesPosition[BestColonyInd, ,drop = FALSE]
  TheEmpire$ImperialistCost <- TheEmpire$ColoniesCost[BestColonyInd]
  TheEmpire$ImperialistInnerParam <- TheEmpire$ColoniesInnerParam[BestColonyInd, ,drop = FALSE]

  TheEmpire$ColoniesPosition[BestColonyInd,] <- OldImperialistPosition
  TheEmpire$ColoniesCost[BestColonyInd] <- OldImperialistCost
  TheEmpire$ColoniesInnerParam[BestColonyInd,] <- OldImperialistInnerParam

  }
  return(TheEmpire)
}
