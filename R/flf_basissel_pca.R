#' Helper function to perform principal component analysis.
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

  Decode <- function(Ystar, k) {
    Yhat_cent <- crossprod(t(Ystar), t(as.matrix(phi_t[, 1:k])))
    Yhat <- sweep(Yhat_cent, MARGIN = 2, STATS = mu_t, FUN = "+")
    Yhat
  }
  return(list(mu_t = mu_t, phi_t = phi_t, Encode = compiler::cmpfun(Encode), Decode = compiler::cmpfun(Decode)))
}


#' Assess losslessness of PCA latent feature representation.
#' @param mat an n times p data matrix.
#' @param kf k, the number of folds in k-fold cross-validation,
#' @param center Whether or not the data should be mean-centered prior to performing the latent feature transformation.
#' @param scale Whether the data should be scaled prior to analysis, so that each variable has standard deviation 1.
#'
#' @return
#' @export
#'
#' @examples
flf_basissel_pca <- function(mat, kf, latent_dim_from = 1, latent_dim_to = min(ncol(mat) - 1, nrow(mat) - 1), latent_dim_by = 1, verbose = TRUE) {

  #* SET UP: MATRIX SIZE AND FUNCTIONS
  mat <- as.matrix(mat)
  n <- nrow(mat)
  p <- ncol(mat)

  breaks <- seq(latent_dim_from, latent_dim_to, by = latent_dim_by)
  breaks <- breaks[which(breaks <= min(n - 1, p - 1))]

  LearnOut <- learn_pca(mat)
  breaks <- breaks[which(breaks <= ncol(LearnOut[["phi_t"]]))]
  q <- length(breaks)

  Encode <- LearnOut[["Encode"]]
  Decode <- LearnOut[["Decode"]]

  corM_t <- RESS <- vector(mode = "numeric", length = q)

  if (verbose) print("====== Training ======")
  for (j in 1:q) { # q is the maximum number of eigenvectors
    if (verbose) print(paste("= Latent Dim. =", breaks[j]))
    proj_t <- Decode(Encode(mat, breaks[j]), breaks[j])
    corM_t[j] <- get_one_minus_squared_correlation(observed = c(mat), predicted = c(proj_t))
  }

  #* CROSS VALIDATION
  if (verbose) print(paste0("====== Performing ", kf, "-fold CV ======="))

  n_cv <- floor(n - n / kf)

  r <- max(which(breaks <= min(n_cv, p)))

  proj_v <- array(NA, c(n, p, q)) # to store the cross-validation predictions
  rho_v <- array(NA, c(n, q)) # to store the individual observations correlations between cross-validation predictions and data for each latent dimension
  Qrho_v <- array(NA, c(n, q)) # version of `rho_v` with it's rows sorted separately for each column (for heatmap)
  PRESS <- corM_v <- vector(mode = "numeric", length = q) # overall/ total correlations/ residual sum of squares between cross-validation predictions and data for each latent dimension
  PRESS[1:q] <- corM_v[1:q] <- NA

  mat <- mat[sample(n), ] # shuffle
  folds <- cut(seq(1, n), breaks = kf, labels = FALSE)

  for (i in 1:kf) {
    # Loop through folds:
    if (verbose) print(paste("==== Fold ", i, "===="))
    kind <- which(folds == i, arr.ind = TRUE)
    mati <- mat[kind, ]
    LearnOutv <- learn_pca(mat[-kind, ])
    Encodev <- LearnOutv[["Encode"]]
    Decodev <- LearnOutv[["Decode"]]
    for (j in 1:r) {
      # Loop through Latent Dimension Sizes:
      if (verbose) print(paste("= Latent Dim. =", breaks[j]))
      proj_v[kind, , j] <- Decodev(Encodev(mati, breaks[j]), breaks[j])
    }
  }

  # Compute Summary Measures on CV Results: ---------------------------------
  if (verbose) print(paste0("====== Finished ", kf, "-fold CV, Summarising Results ======="))

  ## Overall Measures  ------------------------------------------------------
  for (j in seq_len(r)) {
    corM_v[j] <- get_one_minus_squared_correlation(observed = c(mat), predicted = c(proj_v[, , j]))
    PRESS[j] <- norm(as.matrix(proj_v[, , j] - mat), type = "F")^2
  }

  ## Individual-Level Measures ----------------------------------------------
  for (l in seq_len(r)) {
    # loop through latent dimension sizes.
    x <- as.matrix(proj_v[, , l])
    rho_v[, l] <- sapply(
      seq.int(dim(x)[1]),
      function(i) get_one_minus_squared_correlation(observed = as.matrix(mat)[i, ], predicted = as.matrix(x)[i, ])
    )
  }
  rho_v <- matrix(unlist(rho_v), n, q)

  # Sort rows of correlation matrix
  for (s in seq_len(r)) {
    Qrho_v[, s] <- sort(rho_v[, s])
  }


  #* OUTPUT
  out <- list(corM_t = corM_t, corM_v = corM_v, rho_v = rho_v, Qrho_v = Qrho_v, breaks = breaks, r = r, q = q, n = n, p = p, Encode = Encode, Decode = Decode)
}

flf_basissel_pca <- compiler::cmpfun(flf_basissel_pca)
