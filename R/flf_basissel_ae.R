#' Learning function for AE
#'
#' @param Y an n times p data matrix.
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom magrittr "%>%"

learn_ae <- function(Y, k, ae_args) {

  # unpack arguments:
  layer_1_dim <- ifelse(is.null(ae_args[["layer_1_dim"]]), 600, ae_args[["layer_1_dim"]])
  layer_2_dim <- ifelse(is.null(ae_args[["layer_2_dim"]]), 200, ae_args[["layer_2_dim"]])
  link_fun <- ifelse(is.null(ae_args[["link_fun"]]), "sigmoid", ae_args[["link_fun"]])
  epochs <- ifelse(is.null(ae_args[["epochs"]]), 100, ae_args[["epochs"]])
  loss <- ifelse(is.null(ae_args[["loss"]]), "mean_squared_error", ae_args[["loss"]])
  batch_size <-  ifelse(is.null(ae_args[["batch_size"]]), 16, ae_args[["batch_size"]])


  if (!(link_fun %in% c("sigmoid", "linear"))) stop("Link function must be either linear or sigmoid.")
  if (!(loss %in% c("mean_squared_error", "binary_crossentropy"))) stop("Loss must be mean_squared_error or binary_crossentropy")

  p <- ncol(Y)

  # Define the encoder
  encoder <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = layer_1_dim, activation = "relu", input_shape = p) %>%
    keras::layer_dense(units = layer_2_dim, activation = "relu") %>%
    keras::layer_dense(units = k, activation = "linear", name = "bottleneck")

  # Define the decoder
  decoder <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = layer_2_dim, activation = "relu", input_shape = k) %>%
    keras::layer_dense(units = layer_1_dim, activation = "relu") %>%
    keras::layer_dense(units = ncol(Y), activation = link_fun)

  # Connect them to create the autoencoder
  autoencoder <- keras::keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder %>% keras::compile(optimizer = "adam", loss = loss)
  autoencoder %>% keras::fit(Y,
    Y,
    epochs = epochs,
    batch_size = batch_size,
    verbose = FALSE
  )

  Extract <- function(Y) {
    predict(encoder, Y, verbose = FALSE)
  }

  Transform <- function(Ystar) {
    predict(decoder, Ystar, verbose = FALSE)
  }

  return(list(Extract = Extract, Transform = Transform))
}

#' Basis selection for AE.
#'
#' @param Y an n times p data matrix.
#'
#' @return
#' @export
#'
#' @examples
flf_basissel_ae <- function(mat, kf, lim = lim, incr = incr, ae_args = list()) { # check default for breaks

  #* SET UP: MATRIX SIZE AND FUNCTIONS
  n <- nrow(mat)
  p <- ncol(mat)

  #* TRAINING - ALL DATA
  breaks <- c(seq(1, lim, by = incr))
  q <- r <- length(breaks)

  diff_t <- array(0, c(n, q))
  corM_t <- rep(0, q)
  RESS <- rep(0, q)
  mse <- rep(0, q)

  print("====== Training ======")
  for (j in 1:q) {
    print(paste("= Latent Dim. =", j))
    learnout_j <- learn_ae(Y = mat, k = breaks[j], ae_args = ae_args)
    Extract_j <- learnout_j[["Extract"]]
    Transform_j <- learnout_j[["Transform"]]
    proj_t <- Transform_j(Extract_j(mat))
    corM_t[j] <- cor(c(proj_t), c(mat))^2
  }

  #* CROSS VALIDATION
  print(paste0("====== Performing ", kf, "-fold CV ======="))
  proj_v <- array(0, c(n, p))
  # r <- ifelse(n < p, max(floor(q-q/kf), 2), q)
  # r <- ifelse(r<10, r, 10)

  rho_v <- array(0, c(n, r))
  Qrho_v <- array(0, c(n, r))
  diff_v <- array(0, c(n, r))
  corM_v <- rep(0, r)
  PRESS <- c(1:r)
  W <- c(1:(r - 1))

  mat <- mat[sample(n), ]
  folds <- cut(seq(1, n), breaks = kf, labels = FALSE)

  for (i in 1:kf) {
    print(paste("==== Fold ", i, "===="))
    kind <- which(folds == i, arr.ind = TRUE)
    mati <- mat[kind, ]

    for (j in 1:r) {
      print(paste("= Latent Dim. =", j))
      learnout_v_j <- learn_ae(mat[-kind, ], k = breaks[j], ae_args = ae_args)
      Extract <- learnout_v_j[["Extract"]]
      Transform <- learnout_v_j[["Transform"]]

      proj_v[kind, ] <- Transform(Extract(mati))
      proji <- as.matrix(if (NROW(proj_v[kind, ]) > NCOL(proj_v[kind, ]) && (n < p || NROW(proj_v[kind, ]) == p)) t(proj_v[kind, ]) else proj_v[kind, ])
      matir <- if (NROW(mati) < NCOL(mati) || (n > p & NCOL(mati) == p)) as.matrix(t(mati)) else as.matrix(mati)
      rho_v[kind, j] <- sapply(
        seq.int(dim(proji)[1]),
        function(i) cor(as.matrix(proji)[i, ], as.matrix(t(matir))[i, ])^2
      )
      corM_v[j] <- cor(c(proj_v), c(mat))^2
      PRESS[j] <- norm(as.matrix(proj_v - mat), "F")^2
    }
  }

  # sort rows of correlation matrix
  for (s in 1:ncol(rho_v)) {
    Qrho_v[, s] <- sort(rho_v[, s], na.last = TRUE)
  }

  vline <- 0
  out <- list(breaks = breaks, corM_t = corM_t, rho_v = rho_v, Qrho_v = Qrho_v, vline = vline, r = r, n = n, p = p, Extract = Extract, Transform = Transform)
}
