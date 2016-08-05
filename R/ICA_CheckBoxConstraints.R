CheckBoxConstraints <- function(Position, Lower, Upper, Strategy){
  # Position: matrix of positions
  # Lower: vector of lower bound corresponding to each row
  # Upper: upper bound
  # Startegy: character can be 'perturbed' or 'random'
  # output: the corrected positions

  Position <- Position

  n_row <- dim(Position)[1]

  if(Strategy == "perturbed"){
    ExceedLower <- t(sapply(1:n_row, function(j) Position[j,, drop = FALSE] < Lower))
    if(dim(Position)[2] == 1) ##because it gives us wrong answer when we have a matrix whit one colmun
      ExceedLower <- t(ExceedLower)

    if(any(ExceedLower)){
      ExceedLowerIndex <- which(ExceedLower, arr.ind = TRUE)
      L <-  Lower[ExceedLowerIndex[, 2]]
      U <- Upper[ExceedLowerIndex[, 2]]
      x <- Position[ExceedLowerIndex]
      Position[ExceedLowerIndex] <- 2*L - (x + floor((L-x)/(U-L))*(U-L))
      L <- U <- x <- NA
    }


    ExceedUpper <- t(sapply(1:n_row, function(j) Position[j, , drop = FALSE] > Upper))
    if(dim(Position)[2] == 1) ##because it gives us wrong answer when we have a matrix whit one colmun
      ExceedUpper <- t(ExceedUpper)
    if(any(ExceedUpper)) {
      ExceedUpperIndex <- which(ExceedUpper, arr.ind = TRUE)
      L <-  Lower[ExceedUpperIndex[, 2]]
      U <- Upper[ExceedUpperIndex[, 2]]
      x <- Position[ExceedUpperIndex]
      Position[ExceedUpperIndex] <- 2*U - (x - floor((x-U)/(U-L))*(U-L))
    }
  }

  if(Strategy == "random"){

    ExceedLower <- t(sapply(1:n_row, function(j) Position[j, , drop = FALSE] < Lower))
    if(dim(Position)[2] == 1) ##because it gives us wrong answer when we have a matrix whit one colmun
      ExceedLower <- t(ExceedLower)
    if(any(ExceedLower)) {
      ExceedLowerIndex <- which(ExceedLower, arr.ind = TRUE)
      Position[ExceedLowerIndex] <- Lower[ExceedLowerIndex[, 2]] + runif(dim(ExceedLowerIndex)[1])
    }

    ExceedUpper <- t(sapply(1:n_row, function(j) Position[j, , drop = FALSE] > Upper))
    if(dim(Position)[2] == 1) ##because it gives us wrong answer when we have a matrix whit one colmun
      ExceedUpper <- t(ExceedUpper)
    if(any(ExceedUpper)) {
      ExceedUpperIndex <- which(ExceedUpper, arr.ind = TRUE)
      Position[ExceedUpperIndex] <- Upper[ExceedUpperIndex[, 2]] - runif(dim(ExceedUpperIndex)[1])
    }
  }

  return(Position)
}

