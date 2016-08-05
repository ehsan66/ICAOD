#' @importFrom stats runif
GenerateNewCountry <- function(NumOfCountries,
                               lower,
                               upper,
                               sym,
                               sym_point,
                               x_id,
                               w_id,
                               n_independent,
                               equal_weight){

  # sym: if symetric designs should be created
  # sym_point: a point that the design should eb symetrix with respect to
  # n_independent: number of independnt variables. We must know the dimesnion of the design space.
  # when 'sym = TRUE', you need 'x_id' and 'w_id'. 'x_id' is the vector of the column indices that contain the design points
  # for eaxmple if k = 5 then we have, then point 3 is the sym_point and is known. so we should find 4 position for other points and 5 positions for corresponding weights
  # the indices in this case is x_id <- 1:4 and w_id <- 5:9

  # in the current code x_id and w_id is calculated in updtae_ICA as followong
  ##   if(AlgorithmParams$Symetric)
  ##    x_id <- 1:floor(ProblemParams$k/2) else
  ##       x_id <- 1:(ProblemParams$k * AlgorithmParams$n_independent)
  ##    w_id <- (x_id[length(x_id)] + 1):length(ProblemParams$VarMinOuter)

  # output is a matrix of initial solution. each row is points and weights of the design.

  ##each column is one country
  VarMinMatrix <-  matrix(lower, NumOfCountries, length(lower), byrow = TRUE)
  VarMaxMatrix <-  matrix(upper, NumOfCountries, length(upper), byrow = TRUE)
  NewCountry <- (VarMaxMatrix - VarMinMatrix) * matrix(runif (length(VarMaxMatrix)), dim(VarMaxMatrix)[1],  dim(VarMaxMatrix)[2]) + VarMinMatrix

  ## we are sorting even a non symetric design!
  npoint <- length(x_id)


  ### sort only if the number of independent variables is 1, otherwise it's a bug!
  if (n_independent == 1)
    NewCountry[,1:npoint] <- t(apply(NewCountry[,1:npoint, drop = FALSE], 1, sort))
  ##we dont sort the weights based on the points becuase it is an inital random country!
  if(!equal_weight){
    w_mat <- NewCountry[, w_id, drop = FALSE]
  ##ony for symmetric case is useful!
  if (sym)
    even_odd <- ifelse(length(lower)%% 2 == 0, "even", "odd") else
      even_odd <- NA
  ## the sum of weight will be one!!
  w_sum_to_one <- SumToOne(w_mat = w_mat, sym = sym, even_odd = even_odd)
  NewCountry[, w_id] <- w_sum_to_one
  }

  return(NewCountry)
}



