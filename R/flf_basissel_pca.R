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

  svdY <- corpcor::fast.svd(Y) # change this back to fast.svd if we can install corpor (right now not available for R 4.0)
  phi_t <- svdY$v
  Extract <- function(Y, k) {
    Y <- if (NROW(Y) < NCOL(Y) || (n > p & NCOL(Y) == p)) as.matrix(t(Y)) else as.matrix(Y)
    Ystar <- crossprod(as.matrix(Y), as.matrix(phi_t[, 1:k]))
  }
  Transform <- function(Ystar, k) {
    Yhat <- crossprod(t(Ystar), t(as.matrix(phi_t[, 1:k])))
  }
  return(list(phi_t = phi_t, Extract = compiler::cmpfun(Extract), Transform = compiler::cmpfun(Transform)))
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
flf_basissel_pca <- function(mat, kf, lim = min(ncol(mat) - 1, nrow(mat) - 1), incr = 1, center = TRUE, verbose = TRUE) {
  set.seed(123)
  #* SET UP: MATRIX SIZE AND FUNCTIONS
  mat <- as.matrix(mat)
  n <- nrow(mat)
  p <- ncol(mat)

  breaks <- c(seq(1, lim, by = incr))
  breaks <- breaks[which(breaks <= min(n - 1, p - 1))]

  rep.col <- function(x, n) { # this function repeats a column x n times
    matrix(rep(x, each = n), ncol = n, byrow = TRUE)
  }

  #* REMOVE ROW MEANS
  if (center) {
    mat <- as.matrix(mat - rep.col(apply(mat, 1, mean), p))
  }

  #* TRAINING - ALL DATA
  LearnOut <- learn_pca(mat)

  breaks <- breaks[which(breaks <= ncol(LearnOut[["phi_t"]]))]
  q <- length(breaks)


  Extract <- LearnOut[["Extract"]]
  Transform <- LearnOut[["Transform"]]

  diff_t <- array(0, c(n, q))
  corM_t <- rep(0, q)
  RESS <- rep(0, q)

  if(verbose) print("====== Training ======")
  for (j in 1:q) { # q is the maximum number of eigenvectors
      if(verbose) print(paste("= Latent Dim. =", breaks[j]))
      proj_t <- Transform(Extract(mat, breaks[j]), breaks[j])
      corM_t[j] <- cor(c(proj_t), c(mat))^2
    }

  #* CROSS VALIDATION
  if(verbose) print(paste0("====== Performing ", kf, "-fold CV ======="))

  n_cv <- floor(n - n/kf)

  r <- max(which(breaks <= min(n_cv, p)))

  proj_v <- array(0, c(n, p, r))

  rho_v <- array(0, c(n, r))
  Qrho_v <- array(0, c(n, r))
  diff_v <- array(0, c(n, r))
  corM_v <- rep(0, r)
  PRESS <- c(1:r)
  W <- c(1:(r - 1))

  mat <- mat[sample(n), ] # shuffle
  folds <- cut(seq(1, n), breaks = kf, labels = FALSE)

  for (i in 1:kf) {
    if(verbose) print(paste("==== Fold ", i, "===="))
    kind <- which(folds == i, arr.ind = TRUE)
    mati <- mat[kind, ]
    LearnOut <- learn_pca(mat[-kind, ])
    Extractv <- LearnOut[["Extract"]]
    Transformv <- LearnOut[["Transform"]]

    for (j in 1:r) {
      if(verbose) print(paste("= Latent Dim. =", breaks[j]))
      proj_v[kind, , j] <- Transformv(Extractv(mati, breaks[j]), breaks[j])
      proji <- as.matrix(if (NROW(proj_v[kind, , j]) > NCOL(proj_v[kind, , j]) && (n < p || NROW(proj_v[kind, , j]) == p)) {
        t(proj_v[kind, , j])
      } else {
        proj_v[kind, , j]
      })
      matir <- if (NROW(mati) < NCOL(mati) || (n > p & NCOL(mati) == p)) as.matrix(t(mati)) else as.matrix(mati)
      corM_v[j] <- cor(c(proj_v[, , j]), c(mat))^2
      PRESS[j] <- norm(as.matrix(proj_v[, , j] - mat), "F")^2
    }
  }


  for (l in 1:r) {
    x <- as.matrix(proj_v[, , l])
    rho_v[, l] <- sapply(
      seq.int(dim(x)[1]),
      function(i) cor(as.matrix(x)[i, ], as.matrix(mat)[i, ])^2
    )
  }
  rho_v <- matrix(unlist(rho_v), n, r)

  # sort rows of correlation matrix
  for (s in 1:ncol(rho_v)) {
    Qrho_v[, s] <- sort(rho_v[, s])
  }

  #* W STATISTIC
  for (m in 2:r) {
    W[m - 1] <- ((PRESS[m - 1] - PRESS[m]) / (n + p - 2 * m)) / (PRESS[m] / (p * (n - 1) - m * (n + p - m - 1)))
  }
  # vertical line at W statistic
  vline <- ifelse(sum(which(W > 1)) > 0, max(which(W > 1)), 0)


  #* OUTPUT
  out <- list(corM_t = corM_t, rho_v = rho_v, Qrho_v = Qrho_v, vline = vline, breaks = breaks, r = r, q = q, n = n, p = p, Extract = Extract, Transform = Transform)
}

flf_basissel_pca <- compiler::cmpfun(flf_basissel_pca)
