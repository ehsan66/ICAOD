## this files includes the information matrix of the models that have more than one independent variable.
##the C++ information matrix is called here and all of the independent variables converted into one design!!



# Competetive_inhibition
# Noncompetetive_inhibition
# Mixed_type_inhibition
# Uncompetetive_inhibition

FIM_comp_inhibition_x <- function(x, w, param){
  npoint <- length(x)/2
  S <- x[1:npoint]
  I <- x[(npoint+1):(npoint*2)]
  out <- FIM_comp_inhibition(S = S, I = I, w = w, param = param)
  return(out)
}


FIM_noncomp_inhibition_x <- function(x, w, param){
  npoint <- length(x)/2
  S <- x[1:npoint]
  I <- x[(npoint+1):(npoint*2)]
  out <- FIM_noncomp_inhibition(S = S, I = I, w = w, param = param)
  return(out)
}


FIM_mixed_inhibition_x <- function(x, w, param){
  npoint <- length(x)/2
  S <- x[1:npoint]
  I <- x[(npoint+1):(npoint*2)]
  out <- FIM_mixed_inhibition(S = S, I = I, w = w, param = param)
  return(out)
}



FIM_uncomp_inhibition_x <- function(x, w, param){
  npoint <- length(x)/2
  S <- x[1:npoint]
  I <- x[(npoint+1):(npoint*2)]
  out <- FIM_uncomp_inhibition(S = S, I = I, w = w, param = param)
  return(out)
}

