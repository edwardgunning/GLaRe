#' Learning function for principal component analysis (PCA).
#'
#' @param Y an n times p data matrix.
#'
#' @return
#' @export
#'
#' @examples
learn_pca <- function(Y) {
  n <- nrow(Y) # number of observations
  p <- ncol(Y) # number of variables

  mu_t <- apply(Y, 2, mean)
  Y_cent <- sweep(Y, MARGIN = 2, STATS = mu_t, FUN = "-")

  svdY <- corpcor::fast.svd(Y_cent) # change this back to fast.svd if we can install corpor (right now not available for R 4.0)
  phi_t <- svdY$v

  Encode <- function(Y, k) {
    Y_cent <- sweep(Y, MARGIN = 2, STATS = mu_t, FUN = "-")
    Y_cent <- if (NROW(Y_cent) < NCOL(Y_cent) || (n > p & NCOL(Y_cent) == p)) as.matrix(t(Y_cent)) else as.matrix(Y_cent)
    Ystar <- crossprod(as.matrix(Y_cent), as.matrix(phi_t[, 1:k]))
    Ystar
  }

  Decode <- function(Ystar) {
    k <- ncol(Ystar)
    Yhat_cent <- crossprod(t(Ystar), t(as.matrix(phi_t[, 1:k])))
    Yhat <- sweep(Yhat_cent, MARGIN = 2, STATS = mu_t, FUN = "+")
    Yhat
  }
  return(list(mu_t = mu_t, phi_t = phi_t, Encode = compiler::cmpfun(Encode), Decode = compiler::cmpfun(Decode)))
}
