make_vertices <- function(lower, upper){
  #if (length(lower) !=  1){
  ## give the lower and upper of region of uncertainty and the output is the vertices. each row is a vertex
  par_list <- vector("list",length(lower))
  for(i in 1:length(lower))
    par_list[[i]] <- c(lower[i], upper[i])

  vertices <- matrix(unlist(expand.grid(par_list)), nrow = 2^length(upper))
  #} else
    #vertices <-  matrix(c(lower, upper), ncol = 1)
  return(vertices)
}



