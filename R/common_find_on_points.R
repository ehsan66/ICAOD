#... is the further arguments to be passed
##x and w is the arguments that will be passed
find_on_points <- function(fn, ..., points){
  # output calculate fn(points)
  if(!is.function(fn))
    stop("'fn' must be a function.")
 cost <- apply(points, 1, fn,...)

  points <- cbind(points, cost, deparse.level = 0)
##warnings:this is not minima!

  return(list(minima = points, counts = nrow(points)))
}






