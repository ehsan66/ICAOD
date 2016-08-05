ImperialisticCompetition <- function(Empires, Zeta){
  # Empires: all empires
  # perform a imperialistic competition
  # Zeta: a cofficient needed to compute the total cost
  # updating 'ColoniesPosition', "ColoniesCost" and 'ColoniesInnerParam'


  NumOfEmpires <- length(Empires)

  # competition pressure
  if (NumOfEmpires == 1 || runif(1)>.11){
    return(Empires)
  }else{

    #we first select the the empire that possess the weakest colony from weakest empire.
    TotalCosts <- sapply(Empires, "[[", "TotalCost", simplify =TRUE )
    WeakestEmpireInd <- which.max(TotalCosts)

    TotalPowers <- TotalCosts[WeakestEmpireInd] - TotalCosts
    if (all(TotalPowers == 0))
      return(Empires)
    PossessionProbability <- TotalPowers / sum(TotalPowers)

    SelectedEmpireInd <- sample(1:NumOfEmpires, 1, prob = PossessionProbability)


    ## nn is the number of colonies in the weakest Empire. we should choose a colony by random!
    nn <- length(Empires[[WeakestEmpireInd]]$ColoniesCost)
    jj <- sample(1:nn,1)
    ## jj is the index of the colony that should be moved!


    ##add the colony to the selected empire
    Empires[[SelectedEmpireInd]]$ColoniesPosition <- rbind(Empires[[SelectedEmpireInd]]$ColoniesPosition,
                                                           Empires[[WeakestEmpireInd]]$ColoniesPosition[jj, ])

    Empires[[SelectedEmpireInd]]$ColoniesCost <- c(Empires[[SelectedEmpireInd]]$ColoniesCost,
                                                   Empires[[WeakestEmpireInd]]$ColoniesCost[jj])

    Empires[[SelectedEmpireInd]]$ColoniesInnerParam <- rbind(Empires[[SelectedEmpireInd]]$ColoniesInnerParam,
                                                             Empires[[WeakestEmpireInd]]$ColoniesInnerParam[jj, ])

    ##now we should remove the colony from the weakest empire (the jjth colony must be removed).
    Empires[[WeakestEmpireInd]]$ColoniesPosition <- Empires[[WeakestEmpireInd]]$ColoniesPosition[-jj, , drop = FALSE]
    Empires[[WeakestEmpireInd]]$ColoniesCost <- Empires[[WeakestEmpireInd]]$ColoniesCost[-jj]

    #Collapse of the the weakest colony-less Empire
    #if the weakest empire has only one colony we should collapse it and combine it with the strongest empire again
    nn <- length(Empires[[WeakestEmpireInd]]$ColoniesCost)
    if (nn<=1){
      Empires[[SelectedEmpireInd]]$ColoniesPosition <- rbind(Empires[[SelectedEmpireInd]]$ColoniesPosition,
                                                             Empires[[WeakestEmpireInd]]$ImperialistPosition)

      Empires[[SelectedEmpireInd]]$ColoniesCost <- c(Empires[[SelectedEmpireInd]]$ColoniesCost,
                                                     Empires[[WeakestEmpireInd]]$ImperialistCost)

      Empires[[SelectedEmpireInd]]$ColoniesInnerParam <- rbind(Empires[[SelectedEmpireInd]]$ColoniesInnerParam,
                                                             Empires[[WeakestEmpireInd]]$ColoniesInnerParam)
      Empires[[SelectedEmpireInd]]$TotalCost <- Empires[[SelectedEmpireInd]]$ImperialistCost + Zeta * mean(Empires[[SelectedEmpireInd]]$ColoniesCost)
      ##remove the empire from the list (collpase the empire because it has only one colony)
      Empires[[WeakestEmpireInd]] <- NULL
    }
    return(Empires)
  }
}


