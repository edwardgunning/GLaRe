#' Assess losslessness of general feature representation.
#' @param mat an n times p data matrix.
#' @param kf k, the number of folds in k-fold cross-validation,
#' @param center Whether or not the data should be mean-centered prior to performing the latent feature transformation.
#' @param scale Whether the data should be scaled prior to analysis, so that each variable has standard deviation 1.
#'
#' @return
#' @export
#'
#' @examples
flf_basissel <- function(mat, learn, ae_args, kf, latent_dim_from = 1, latent_dim_to = min(ncol(mat), nrow(mat) - 1), latent_dim_by = 1, loss_function = get_one_minus_squared_correlation, learn_function = NULL, verbose = TRUE) {

  # add arguments check(s).
  if(learn == "user" & is.null(learn_function)) stop("learn_function must be supplied if learn = 'user'.")
  if(learn == "user") {
    if(class(learn_function) != "function") stop("learn_function() must be a function.")
  }

  # Data matrix and its dimensions. -----------------------------------------
  if(learn == "dwt.2d") {
    if (!(class(mat) == "array" & length(dim(mat)) == 3)) {
      stop("When learn = 'dwt.2d', mat must be a 3-dimensional array, where each slice is an image.")
    }
    arr <- mat
    arr_dims <- dim(arr)
    n <- arr_dims[1]
    p1 <- arr_dims[2]
    p2 <- arr_dims[3]
    p <- p1 * p2
    if(latent_dim_to == p1) latent_dim_to <- p # temporary fix for now.
    mat <- matrix(arr, nrow = n, ncol = p)
  } else {
    if(!("matrix" %in% class(mat))) stop("mat must be an n times p data matrix")
    n <- nrow(mat)
    p <- ncol(mat)
  }

  # Set up sequence for latent dimensions: ----------------------------------
  breaks <- seq(from = latent_dim_from, to = latent_dim_to, by = latent_dim_by)
  breaks <- if(learn == "pca") {
    breaks[which(breaks <= min(n - 1, p))]
  } else {
    breaks[which(breaks <= p)]
  }

  # Set up feature learning function: ---------------------------------------
  if (learn == "pca") {
    learn_function <- learn_pca
  } else if (learn == "dwt") {
    learn_function <- learn_dwt
  } else if (learn == "dwt.2d") {
    learn_function <- function(Y) {
      learn_dwt.2d(Y = Y, p1 = p1, p2 = p2)
      }
    } else if (learn == "ae") {
    learn_function <- function(Y, k) {
      learn_ae(Y = Y, k = k, ae_args = ae_args)
    }
  }

  # Do training: ------------------------------------------------------------
  if (learn %in% c("pca", "dwt", "dwt.2d")) {
    LearnOut <- learn_function(mat)
    Encode <- LearnOut[["Encode"]]
    Decode <- LearnOut[["Decode"]]
  }

  if(learn == "pca") {
    if(!all(breaks <= ncol(LearnOut[["phi_t"]]))) print("Number of non-zero eigenvectors is less than latent_dim_to")
    breaks <- breaks[which(breaks <= ncol(LearnOut[["phi_t"]]))]
  }

  q <- length(breaks)

  # Set up vectors to score training results: -------------------------------
  corM_t <- vector(mode = "numeric", length = q)

  # Do training: ------------------------------------------------------------
  if (verbose) print("====== Training ======")

  for (j in 1:q) {
    if (verbose) print(paste("= Latent Dim. =", breaks[j]))
    if (learn %in% c("pca", "dwt", "dwt.2d")) { # if learning is agnostic to latent dimension
      proj_t <- Decode(Encode(Y = mat, k = breaks[j]))
    } else { # if learning depends on latent dimension
      learnout_j <- learn_function(Y = mat, k = breaks[j])
      proj_t <- learnout_j[["Decode"]](learnout_j[["Encode"]](mat))
    }
    corM_t[j] <- loss_function(observed = c(mat), predicted = c(proj_t))
    if(learn == "ae") {
      keras::k_clear_session()
      rm(learnout_j)
      gc(verbose = FALSE)
    }
    rm(proj_t)
    gc(verbose = FALSE)
  }





  # Do Cross-Validation: ----------------------------------------------------
  if (verbose) print(paste0("====== Performing ", kf, "-fold CV ======="))

  n_cv <- floor(n - n / kf) # sample sizes for each cv training

  # dimensionality for cross validation.
  r <- if(learn == "pca") {
    max(which(breaks <= min(n_cv - 1, p))) # maximum dimensionality for PCA.
  } else {
    max(which(breaks <= p))
    }

  proj_v <- array(NA, c(n, p, q)) # to store the cross-validation predictions
  rho_v <- array(NA, c(n, q)) # to store the individual observations correlations between cross-validation predictions and data for each latent dimension
  Qrho_v <- array(NA, c(n, q)) # version of `rho_v` with it's rows sorted separately for each column (for heatmap)
  PRESS <- corM_v <- vector(mode = "numeric", length = q) # overall/ total correlations/ residual sum of squares between cross-validation predictions and data for each latent dimension
  PRESS[1:q] <- corM_v[1:q] <- NA

  # shuffle data prior to cross-validation.
  shuffle_inds <- sample(n)
  mat <- mat[shuffle_inds, ]
  folds <- cut(seq_len(n), breaks = kf, labels = FALSE)

  for (i in 1:kf) {
    # Loop through folds:
    if (verbose) print(paste("==== Fold ", i, "===="))
    kind <- which(folds == i, arr.ind = TRUE)
    mati <- mat[kind,, drop = FALSE]
    if (learn %in% c("pca", "dwt", "dwt.2d")) {
      LearnOutv <- learn_function(Y = mat[-kind,, drop = FALSE])
      Encodev <- LearnOutv[["Encode"]]
      Decodev <- LearnOutv[["Decode"]]
    }
    for (j in 1:r) {
      # Loop through Latent Dimension Sizes:
      if (verbose) print(paste("= Latent Dim. =", breaks[j]))
      if (learn %in% c("pca", "dwt", "dwt.2d")) {
        proj_v[kind, , j] <- Decodev(Encodev(Y = mati, k = breaks[j]))
      } else {
        learnout_v_j <- learn_function(Y = mat[-kind,, drop = FALSE], k = breaks[j])
        proj_v[kind, , j] <- learnout_v_j[["Decode"]](learnout_v_j[["Encode"]](mati))
      }
      if(learn == "ae") {
        keras::k_clear_session()
        rm(learnout_v_j)
        gc(verbose = FALSE)
        }
    }
  }


  # Compute Summary Measures on CV Results: ---------------------------------
  if (verbose) print(paste0("====== Finished ", kf, "-fold CV, Summarising Results ======="))

  ## Overall Measures  ------------------------------------------------------
  for (j in seq_len(r)) {
    corM_v[j] <- loss_function(observed = c(mat), predicted = c(proj_v[, , j]))
    PRESS[j] <- norm(as.matrix(proj_v[, , j] - mat), type = "F")^2
  }

  ## Individual-Level Measures ----------------------------------------------
  for (l in seq_len(r)) {
    # loop through latent dimension sizes.
    x <- as.matrix(proj_v[, , l])
    rho_v[, l] <- sapply(
      seq.int(dim(x)[1]),
      function(i) loss_function(observed = as.matrix(mat)[i, ], predicted = as.matrix(x)[i, ])
    )
  }
  rho_v <- matrix(unlist(rho_v), n, q)

  # un-shuffle:
  rho_v[shuffle_inds, ] <- rho_v

  # Sort rows of correlation matrix
  for (s in seq_len(r)) {
    Qrho_v[, s] <- sort(rho_v[, s])
  }

  # Return output in list: --------------------------------------------------
  list(corM_t = corM_t, corM_v = corM_v, rho_v = rho_v, Qrho_v = Qrho_v, breaks = breaks, n = n, p = p, r = r, q = q, n = n, p = p, learn = learn)
}
