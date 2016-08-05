#' @importFrom stats optim

# like optim, but with n.seg. divide the uppper-lower into n.seg intervals and start to find the local optima
# with respect to each endpoints.
# The outputs is a list of minima (cost values is the last column) and the number of function evaluations.
optim2 <- function(fn, ..., lower , upper, control = list(factr = sqrt(.Machine$double.eps)), n.seg){
  if(!is.function(fn))
    stop("'fn' must be a function.")
  fn1 <- function(arg)fn(arg, ...)

  u <- (upper - lower)/n.seg

  partition <- t(sapply(X = (0):(n.seg),
                        FUN = function(i) lower + i * u ))

  if(length(lower) == 1)
    partition <- t(partition)
  p0 <- do.call(`expand.grid`,as.data.frame(partition))

  result <- sapply(X = 1:dim(p0)[1],
                   FUN = function(i)optim(par = p0[i, ], fn = fn1,
                                          lower = lower, upper = upper,
                                          method = "L-BFGS-B",
                                         control = list(factr = control$factr)))

  minima <- result[1:2, 1:dim(p0)[1]]

  minima <- matrix(unlist(minima, use.names=FALSE), ncol = length(lower) + 1, byrow = TRUE)
  #minima <- unique(minima, MARGIN = 1)
  minima <- minima[!duplicated(round(minima[, -dim(minima)[2]], 3)), , drop = FALSE]
  counts <- sum(sapply(result[3,], "[[", 1))

  return(list(minima = minima, counts = counts))
}










