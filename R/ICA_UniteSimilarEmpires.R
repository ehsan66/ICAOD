UniteSimilarEmpires <- function(Empires, Zeta, UnitingThreshold, SearchSpaceSize){
  # Unite the similar empires
  # it can be removed



  TheresholdDistance <- UnitingThreshold * base::norm(as.matrix(SearchSpaceSize), type = '2')
  NumOfEmpires <- length(Empires)


  ##because the length(Empire) = NumOfEmpires is dynamic (it will be decreased if uniting happens we should use while)
  counter1 <- 1
  counter2 <- 2

  while(counter1 <= (NumOfEmpires-1)){
    while(counter2 <= NumOfEmpires){
      DistanceVector <- Empires[[counter1]]$ImperialistPosition - Empires[[counter2]]$ImperialistPosition
      Distance <- base::norm(DistanceVector, "2")
      #cat(Distance <= TheresholdDistance, "\n")
      if(Distance <= TheresholdDistance){
        if (Empires[[counter1]]$ImperialistCost < Empires[[counter2]]$ImperialistCost){
          BetterEmpireInd <- counter1
          WorseEmpireInd <- counter2} else{
            BetterEmpireInd <- counter2
            WorseEmpireInd <- counter1
          }

        Empires[[BetterEmpireInd]]$ColoniesPosition <- rbind(Empires[[BetterEmpireInd]]$ColoniesPosition,
                                                        Empires[[WorseEmpireInd]]$ImperialistPosition,
                                                        Empires[[WorseEmpireInd]]$ColoniesPosition)

        Empires[[BetterEmpireInd]]$ColoniesCost <- c(Empires[[BetterEmpireInd]]$ColoniesCost,
                                                    Empires[[WorseEmpireInd]]$ImperialistCost,
                                                    Empires[[WorseEmpireInd]]$ColoniesCost)

        Empires[[BetterEmpireInd]]$ColoniesInnerParam <- rbind(Empires[[BetterEmpireInd]]$ColoniesInnerParam,
                                                               Empires[[WorseEmpireInd]]$ImperialistInnerParam,
                                                               Empires[[WorseEmpireInd]]$ColoniesInnerParam)

        # Update TotalCost for new United Empire
        Empires[[BetterEmpireInd]]$TotalCost = Empires[[BetterEmpireInd]]$ImperialistCost + Zeta * mean(Empires[[BetterEmpireInd]]$ColoniesCost)
        Empires[[WorseEmpireInd]] <- NULL
        NumOfEmpires <- NumOfEmpires - 1
      }
      counter2 <- counter2 + 1
    }
    counter1 <- counter1 + 1
  }




  return(Empires)
}






