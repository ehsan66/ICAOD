# ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
#      dnx <- dimnames(X)
#      if (is.null(dnx))
#          dnx <- vector("list", 2)
#      s <- svd(X)
#      nz <- s$d > tol * s$d[1]
#      structure(if (any(nz))
#          s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
#          else X, dimnames = dnx[2:1])
# }

#### Moore-Penrose Matrix Inverse
mpginv <- function (a, tol = sqrt(.Machine$double.eps))
{

  if (length(dim(a)) > 2 || !(is.numeric(a) || is.complex(a)))
    stop("a must be a numeric or complex matrix")
  if (!is.matrix(a))
    a <- as.matrix(a)
  asvd <- La.svd(a)
  if (is.complex(a)) {
    asvd$u <- Conj(asvd$u)
    asvd$v <- t(Conj(asvd$vt))
  }
  else {
    asvd$v <- t(asvd$vt)
  }
  Positive <- asvd$d > max(tol * asvd$d[1], 0)
  if (!any(Positive))
    array(0, dim(a)[2:1])
  else asvd$v[, Positive] %*% ((1/asvd$d[Positive]) * t(asvd$u[ , Positive]))
}
