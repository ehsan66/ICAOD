

blatex <- function(obj){

  if (obj$arg$prior == "uniform"){
    Theta <- paste("$", paste("[", paste(obj$arg$uniform$min, obj$arg$uniform$max, sep = ", "), "]", sep = "", collapse = "\\times"), "$", sep = "")
  }
  iter <- length(obj$evol)
  crt <- paste("$", round(obj$evol[[iter]]$min_cost, 5), "$", sep = "")
  ELB <- paste("$", round(obj$evol[[iter]]$ELB, 5), "$", sep = "")
  temp <- paste("$", paste(round(obj$evol[[iter]]$x, 5), " (", round(obj$evol[[iter]]$w, 5), ")", sep = ""), "$", sep = "", collapse = "\\\\ ")
  xi <- paste("\\specialcell{" ,temp, "}", sep = "")
  cpu <- "NA"
  #relative <- "NA"
  row <- paste(Theta, xi, crt, ELB, iter, cpu, sep = " & ")
  cat(row)


}

