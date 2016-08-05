CreateInitialEmpires <- function(sorted_Countries,
                                 sorted_Cost,
                                 sorted_InnerParam,
                                 Zeta,
                                 NumOfInitialImperialists,
                                 NumOfAllColonies
                                 ){
  # create initial empires, as a list

  # sorted_Countries: matrix of countries, must be sorted while the 'nimp' fisrt of them will be choosen as imperialists
  #  sorted_Cost: corresponding cost of sorted_count
  # sorted_InnerParam: corresponding inner parameters for each countries
  # Zeta
  # NumOfInitialImperialists: number of imperialists
  # ncolon: number of colonies
  # output: a list contains positions and costs for the corresponding colonies and imperialists.
  AllImperialistsPosition <- sorted_Countries[1:NumOfInitialImperialists, , drop = FALSE]
  AllImperialistsCost = sorted_Cost[1:NumOfInitialImperialists]
  AllImperialistsInnerParam <- sorted_InnerParam[1:NumOfInitialImperialists, , drop = FALSE]


  ## other colonies, the imperialists are excluded
  AllColoniesPosition = sorted_Countries[-(1:NumOfInitialImperialists), , drop = FALSE]
  AllColoniesCost = sorted_Cost[-c(1:NumOfInitialImperialists)]
  AllColoniesInnerParam= sorted_InnerParam[-(1:NumOfInitialImperialists), , drop = FALSE]


  # calculate imperialists power
  if (max(AllImperialistsCost)>0)
    AllImperialistsPower <- 1.3 * max(AllImperialistsCost) - AllImperialistsCost else
      AllImperialistsPower <- 0.7 * max(AllImperialistsCost) - AllImperialistsCost
  AllImperialistsPower <- AllImperialistsPower/sum(AllImperialistsPower)

  # finding the number of colonies for each empire
  ## number of colonies of each impire. each one should at leats has one colony
  AllImperialistNumOfColonies <- rep(1, length(AllImperialistsPower))
  NumOfRemainColonies <- NumOfAllColonies - sum(AllImperialistNumOfColonies)
  AllImperialistNumOfColonies <- floor(AllImperialistsPower * NumOfRemainColonies)  + AllImperialistNumOfColonies
  ##if still some colonies remain add it to the strongest one
  diff_col <- NumOfAllColonies - sum(AllImperialistNumOfColonies)
  # we add or to the strongest imperialist
  AllImperialistNumOfColonies[1] <-  AllImperialistNumOfColonies[1] + diff_col
  ## modfy the number of colonies for last Impire to be sure that the total number of colonies not less or more than ncolon
  ## warning: bug prone, becasue it may delet all the colonies when it only has one
  last_empire_id <- length(AllImperialistNumOfColonies)
  AllImperialistNumOfColonies[last_empire_id] <- NumOfAllColonies - sum(AllImperialistNumOfColonies[-last_empire_id])



  RandomIndex <- sample(x = NumOfAllColonies, size = NumOfAllColonies, replace = FALSE)
  RandonIndexList <- vector("list", NumOfInitialImperialists)
  InitialInd <- 1
  LastInd <- AllImperialistNumOfColonies[1]

  for(jj in 1:NumOfInitialImperialists){
    if( InitialInd > length(RandomIndex))
      RandonIndexList [[jj]] <- 0 else
        RandonIndexList [[jj]] <- RandomIndex[InitialInd:LastInd]
    InitialInd  <- 1 + LastInd
    LastInd <- LastInd +   AllImperialistNumOfColonies[jj + 1]
  }

  ##now we create the empires. we sent the Imperialist  position and cost and to a list
  Empires <- vector("list", NumOfInitialImperialists)
  for(i in 1:length(AllImperialistNumOfColonies)){
    Empires[[i]]$ImperialistPosition <- AllImperialistsPosition[i, , drop = FALSE]
    Empires[[i]]$ImperialistCost <- AllImperialistsCost[i]
    Empires[[i]]$ImperialistInnerParam <- AllImperialistsInnerParam[i, , drop = FALSE]
    if(all(RandonIndexList[[i]]== 0)){
      stop("One initial empire has no colony.\nIncrease the number of colonies or decrease the number of imperialists.\nUsually the number of imperialists set to be 10 percent of the colonies.")
    }else{
      Empires[[i]]$ColoniesPosition <- AllColoniesPosition[RandonIndexList[[i]], , drop = FALSE]
      Empires[[i]]$ColoniesCost <- AllColoniesCost[RandonIndexList[[i]]]
      Empires[[i]]$ColoniesInnerParam <- AllColoniesInnerParam[RandonIndexList[[i]], , drop = FALSE]
      Empires[[i]]$TotalCost <- Empires[[i]]$ImperialistCost + Zeta * mean(Empires[[i]]$ColoniesCost)
    }

  }

  return(Empires)

}

