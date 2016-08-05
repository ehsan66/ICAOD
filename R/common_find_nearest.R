find_nearest <- function(values, tol, compare_with_minimum){
  # return the indices of "values" that are near to minimum or maximim with 'tol'

  if(compare_with_minimum)
    global_index <- which.min(values) else
      global_index <- which.max(values)
  global <- values[global_index]
  keep <- which(values - global < tol)
  return(keep)
}



