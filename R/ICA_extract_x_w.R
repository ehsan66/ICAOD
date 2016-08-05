ICA_extract_x_w <- function(x, w, sym, sym_point){
  ## when the deign is symetrix, extratc the x and w
  include_symetric <- length(x) != length(w)
  x <- c(x, rev(2 * sym_point - x))
  if(include_symetric)
    x <- c(x, sym_point)
  if(include_symetric){
    w_sym <- w[length(w)]
    w <- w[-length(w)]
    w <- c(w, rev(w))
    w <- c(w, w_sym)
  }else
    w <- c(w, rev(w))

  return(list(x=x, w=w))
}


