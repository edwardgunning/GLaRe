#' Assess losslessness of 2-D DWT latent feature representation.
#' @param mat an n times p times q data array, containing p times q images in the n slices
#' @param kf k, the number of folds in k-fold cross-validation,
#' @param center Whether or not the data should be mean-centered prior to performing the latent feature transformation.
#' @param scale Whether the data should be scaled prior to analysis, so that each variable has standard deviation 1.
#'
#' @return
#' @export
#'
#' @examples
flf_basissel_dwt.2d <- function(mat, kf, latent_dim_from = 1, latent_dim_to = min(ncol(mat) - 1, nrow(mat) - 1), latent_dim_by = 1, verbose = TRUE) {
  if (!length(dim(mat)) == 3) {
    stop("mat must be a 3-dimensional array, where each slice is an image.")
  }
  mat_dims <- dim(mat)
  n <- mat_dims[1]
  p1 <- mat_dims[2]
  p2 <- mat_dims[3]
  p <- p1 * p2
  mat_flattened <- keras::array_reshape(mat, dim = c(n, p))

  breaks <- seq(1, latent_dim_to, by = latent_dim_by)
  breaks <- breaks[which(breaks <= min(n - 1, p - 1))]
  q <- length(breaks)

  rep.col <- function(x, n) { # this function repeats a column x n times
    matrix(rep(x, each = n), ncol = n, byrow = TRUE)
  }

  #* TRAINING - ALL DATA
  LearnOut <- learn_dwt.2d(Y = mat)

  Encode <- LearnOut[["Encode"]]
  Decode <- LearnOut[["Decode"]]

  diff_t <- array(0, c(n, q))
  corM_t <- rep(0, q)
  RESS <- rep(0, q)

  if (verbose) print("====== Training ======")
  for (j in 1:q) { # q is the maximum number of eigenvectors
    if (verbose) print(paste("= Latent Dim. =", breaks[j]))
    proj_t <- keras::array_reshape(x = Decode(Ystar = Encode(Y = mat, k = breaks[j])), c(n, p))
    corM_t[j] <- 1 - cor(c(proj_t), c(mat))^2
  }


  #* CROSS VALIDATION
  if (verbose) print(paste0("====== Performing ", kf, "-fold CV ======="))

  n_cv <- floor(n - n / kf)

  r <- max(which(breaks <= min(n_cv, p)))

  proj_v <- array(0, c(n, p, r))

  rho_v <- array(0, c(n, r))
  Qrho_v <- array(0, c(n, r))
  diff_v <- array(0, c(n, r))
  corM_v <- rep(0, r)
  PRESS <- c(1:r)
  W <- c(1:(r - 1))

  mat <- mat[sample(n), , ] # shuffle
  folds <- cut(seq(1, n), breaks = kf, labels = FALSE)

  for (i in 1:kf) {
    if (verbose) print(paste("==== Fold ", i, "===="))
    kind <- which(folds == i, arr.ind = TRUE)
    mati <- mat[kind, , ]
    mati_flattened <- keras::array_reshape(mati, dim = c(length(kind), p))
    LearnOut <- learn_dwt.2d(mat[-kind, , ])
    Encodev <- LearnOut[["Encode"]]
    Decodev <- LearnOut[["Decode"]]

    for (j in 1:r) {
      if (verbose) print(paste("= Latent Dim. =", breaks[j]))
      proj_v[kind, , j] <- keras::array_reshape(Decodev(Ystar = Encodev(Y = mati, k = breaks[j])), dim = c(length(kind), p))
      proji <- as.matrix(if (NROW(proj_v[kind, , j]) > NCOL(proj_v[kind, , j]) && (n < p || NROW(proj_v[kind, , j]) == p)) {
        t(proj_v[kind, , j])
      } else {
        proj_v[kind, , j]
      })
    }
  }


  for (j in 1:r) {
    corM_v[j] <- 1 - cor(c(proj_v[, , j]), c(mat))^2
    PRESS[j] <- norm(as.matrix(proj_v[, , j] - mat_flattened), "F")^2
  }


  for (l in 1:r) {
    x <- as.matrix(proj_v[, , l])
    rho_v[, l] <- sapply(
      seq.int(dim(x)[1]),
      function(i) 1 - cor(as.matrix(x)[i, ], as.matrix(mat_flattened)[i, ])^2
    )
  }
  rho_v <- matrix(unlist(rho_v), n, r)

  # sort rows of correlation matrix
  for (s in 1:ncol(rho_v)) {
    if (all(is.na(rho_v[, s]))) Qrho_v[, s] <- rho_v[, s] else Qrho_v[, s] <- sort(rho_v[, s])
  }


  vline <- NA

  #* OUTPUT
  out <- list(corM_t = corM_t, rho_v = rho_v, Qrho_v = Qrho_v, vline = vline, breaks = breaks, r = r, q = q, n = n, p = p, Encode = Encode, Decode = Decode)
}

flf_basissel_dwt.2d <- compiler::cmpfun(flf_basissel_dwt.2d)
