#' Perform PCA with Encoding and Decoding Functions
#'
#' This function implements Principal Component Analysis (PCA), returning the
#' mean-centered data, principal components, and encoding/decoding functions for
#' dimensionality reduction and reconstruction.
#'
#' @param Y A numeric matrix with `n` rows (observations) and `p` columns (variables).
#'
#' @return A list containing:
#' \itemize{
#'   \item `mu_t`: A vector of column means used for centering the data.
#'   \item `phi_t`: A matrix of eigenvectors (principal components).
#'   \item `Encode`: A function to encode data into the reduced PCA space.
#'   \item `Decode`: A function to reconstruct data from the PCA space.
#' }
#' @examples
#' # Example usage
#' library(GLaRe)
#' Y <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' pca_model <- learn_pca(Y)
#' encoded_data <- pca_model$Encode(Y, k = 3)
#' reconstructed_data <- pca_model$Decode(encoded_data)
#'
#' @export
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
