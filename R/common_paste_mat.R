# cat matrix
paste_mat <- function(mat){
  temp <- ""
  for(i in 1:dim(mat)[1])
    temp <- paste(temp, "\n", paste(mat[i, ], sep = "", collapse = "  "),  sep = "")

  #cat(temp)
  return(temp)
}


