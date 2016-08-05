SumToOne <- function(w_mat, sym, even_odd){
  #even_odd is NA when sym == FALSE
  # even_odd: is the number of support points is odd or even. needed when sym = TRUE
  # because when sym = TRUE the sum opf weights must be one, but symetric

  if(!sym)
    w_mat_out <- w_mat/apply(w_mat, 1, sum) else{
      n_col <- dim(w_mat)[2]

      if(even_odd == "even"){
        all_w <- cbind(w_mat, w_mat[, rev(1:n_col), drop = FALSE])
        all_w <- all_w/apply( all_w, 1, sum)
        w_mat_out <- all_w[, 1:n_col, drop = FALSE]
      }
      if(even_odd == "odd"){
        sym_w <- w_mat[,n_col, drop = FALSE]
        ##now we remove the symetric weight!
        w_mat <- w_mat[, -n_col, drop = FALSE]
        all_w <- cbind(w_mat, w_mat[, rev(1:dim(w_mat)[2]), drop = FALSE])
        all_w <- cbind(all_w, sym_w)
        all_w <- all_w/apply( all_w, 1, sum)
        w_mat_out <- cbind(all_w[,1:(dim(all_w)[2]/2), drop =FALSE], all_w[, dim(all_w)[2], drop = FALSE])
      }
    }



  return(w_mat_out)
}
