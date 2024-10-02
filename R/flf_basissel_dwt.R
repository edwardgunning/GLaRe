# set.seed(123)
# data("glaucoma_images", package = "complex.fundata"); dim(glaucoma.images)
# glauc.data <- t(matrix(glaucoma.images, 120*120, 306))
# glspl <- glauc.data[sample(1:nrow(glauc.data), 20, replace=F), ]



# library(wavelets)
#
# wtData <- NULL
#
# for (i in 1:nrow(mnist_spl)) {
#   a <- as.numeric(mnist_spl[i,])
#   wt <- dwt(a, filter='d4', n.levels=5, boundary='periodic')
#   wtData <- rbind(wtData, unlist(c(wt@W,wt@V[[wt@level]])))
#   }
# wtData <- as.data.frame(wtData)

# https://www.rdocumentation.org/packages/wavelets/versions/0.3-0.2/topics/wt.filter

#*** GET COMPRESS FUNCTION ***#
get_Kj_compress <- function(wavespecs) {
  keep <- wavespecs$keep
  J <- length(wavespecs$Kj)
  Kj_all <- wavespecs$Kj
  Kj <- rep(0, J)

  temp_first <- c(1, 1 + cumsum(Kj_all[1:(J - 1)]))
  temp_last <- cumsum(Kj_all)
  for (i in 1:J) {
    Kj[i] <- length(wavespecs$keep[temp_first[i]:temp_last[i]])
  }
  return(list(Kj, Kj_all))
}

#*** JOINT COMPRESSION - GET EN ***#
Energy <- function(D) {
  En <- D
  total_Energy <- rep(0, nrow(D))
  for (i in 1:nrow(D)) {
    ix <- order(-abs(D[i, ]))
    Csort <- -D[i, ix]
    total_Energy[i] <- sum(Csort^2)
    energy <- cumsum(as.numeric(Csort^2)) / total_Energy[i]
    En[i, ix] <- energy
  }
  return(En)
}


#*** SPLIT MATRIX ***#
split_matrix <- function(M, list_of_rows, list_of_cols) {
  temp <- cbind(list_of_rows, list_of_cols)
  col1 <- 1
  col2 <- list_of_cols[1]
  S <- list()
  for (i in 1:(nrow(temp) - 1)) {
    S[[i]] <- M[1:temp[i, 1], col1:col2]
    col1 <- col1 + temp[i, 2]
    col2 <- col2 + temp[i + 1, 2]
  }
  return(S)
}

#' learn
#' @param
#'
#' @return
#' @export
#'
#' @examples
#*** FEATURE LEARNING FUNCTION ***#
learn_dwt <- function(
    Y,
    wavespecs = list("d4", n.levels = 5, boundary = "periodic")
    # , wavespecs= list(wavelet = "db3", nlevels = c(5,5), compress
    #                   = 1, ndim = 2, t = c(120, 120), rectangular = 1, boundary = c("per", "sym"))
    , breaks = c(seq(0.8, 0.875, by = 0.025), seq(0.9, .95, by = 0.01), seq(.95, 1, by = 0.0004))
    # , breaks
    ) {
  # D <- dwt(as.matrix(Y), wavespecs = ws)
  # wavespecs <- attr(D, "wavespecs")
  wave <- wavelets::dwt(Y, wavespecs[[1]], wavespecs[[2]], wavespecs[[3]]) # , filter = 'd4', n.levels=5, boundary='periodic')
  coeffs <- wave@W
  coeffs[[wave@level]] <- c(wave@W[[wave@level]], wave@V[[wave@level]])
  d <- lapply(coeffs, function(x) matrix(x, nrow(Y), length(x) / nrow(Y)))
  sizes <- lapply(d, function(x) dim(x)[2])
  D <- do.call(cbind, d)

  Dcomp <- Energy(D)
  scree <- apply(Dcomp, 2, mean)
  screeix <- order(scree)
  screeord <- scree[screeix]
  plot(
    screeord,
    type = "l",
    col = "blue",
    lwd = 2,
    xlab = "Wavelet Coefficients",
    ylab = "Scree",
    main = "Mean Energy (Scree)"
  )

  # Extract <- function(Y, k, keep = NULL) {
  #   if (all(is.null(keep))) {
  #     keep <- scree <= breaks[k]
  #   } else {
  #     keep <- keep
  #   }
    Extract <- function(Y, k) {

    keep <- scree <= breaks[k]
    index <- max(which(scree <= breaks[k]))
    D0 <- t(apply(D, 1, function(x) x * keep))

    wave <- wavelets::dwt(Y, wavespecs[[1]], wavespecs[[2]], wavespecs[[3]])
    D1 <- split_matrix(D0, rep(nrow(D), length(sizes) + 1), c(as.numeric(unlist(sizes)), 0)) # $$$$$$$
    D1 <- lapply(D1, function(x) as.matrix(c(x)))
    coeffs <- D1
    coeffs[[wave@level]] <- as.matrix(D1[[wave@level]][1:(length(D1[[length(sizes)]]) / 2)])
    wave@W <- coeffs
    wave@V[[wave@level]] <- as.matrix(D1[[wave@level]][((length(D1[[length(sizes)]]) / 2) + 1):(length(D1[[length(sizes)]]))])
    Ystar <- wave
    return(list(Ystar, keep, index))
  }
  Transform <- function(Ystar) {
    Yhat <- wavelets::idwt(Ystar, fast = TRUE)
  }

  return(list(D = D, Extract = Extract, Transform = Transform))

  # index <- max(which(scree <= breaks[k]))
}


#' DWT
#' @param
#'
#' @return
#' @export
#'
#' @examples
flf_basissel_dwt <- function(
    mat,
    kf,
    learn = learn_dwt,
    wavespecs = list(filter = "d4", n.levels = 5, boundary = "periodic"),
    center = TRUE,
    breaks = c(seq(0.8, 0.975, by = 0.025), seq(0.975, 1, by = 0.001)))
# , breaks=c(seq(0.8, 0.9, by=0.025), seq(0.95, .975, by=0.001), seq(.975, 1, by=0.0001)))
{ # check default for breaks

  #* SET UP: MATRIX SIZE AND FUNCTIONS
  n <- nrow(mat)
  p <- ncol(mat)
  mat_cent <- as.matrix(mat) # no centering

  #* TRAINING - ALL DATA
  learnout <- learn(mat_cent, wavespecs, breaks)
  D <- learnout[[1]]
  # numix <- ncol(D)
  # q <- ifelse(numix<10, numix, 10)
  # breaks <- c(seq(0.5, 1, by = 0.05))
  # breaks <- c(seq(0.5, 0.8, by=0.1), 0.85, seq(0.9, 0.99, by=0.01), 0.999, 1)
  # breaks <- c(seq(0.5, 0.875, by=0.025), seq(0.9, .95, by=0.01), seq(.95, 1, by=0.0004) )
  # breaks <- c(seq(0.8, 0.875, by=0.025), seq(0.9, .95, by=0.01), seq(.95, 1, by=0.0004))
  Extract <- learnout[[2]]
  Transform <- learnout[[3]]
  q <- r <- length(breaks)

  diff_t <- array(0, c(n, q))
  corM_t <- rep(0, q)
  RESS <- rep(0, q)
  indices <- c()

  print("====== Training ======")
  for (j in 1:q) {
    print(paste("= Latent Dim. =", j))
    proj_t <- Transform(Extract(mat, j)[[1]])
    indices <- c(indices, Extract(mat, j)[[2]])
    corM_t[j] <- cor(c(proj_t), c(mat_cent))^2
    RESS[j] <- norm(as.matrix(proj_t - mat_cent), "F")^2
  }
  plot(corM_t, type = "l", col = "blue", lwd = 3)

  #* CROSS VALIDATION
  print(paste0("====== Performing ", kf, "-fold CV ======="))
  proj_v <- array(0, c(n, p))
  # r <- ifelse(n < p, max(floor(q-q/kf), 2), q)
  # r <- ifelse(r<200, r, 200)

  rho_v <- array(0, c(n, r))
  Qrho_v <- array(0, c(n, r))
  diff_v <- array(0, c(n, r))
  corM_v <- rep(0, r)
  PRESS <- c(1:r)
  W <- c(1:(r - 1))

  mat_cent <- mat_cent[sample(n), ]
  folds <- cut(seq(1, n), breaks = kf, labels = FALSE)

  for (i in 1:kf) {
    print(paste("==== Fold ", i, "===="))
    kind <- which(folds == i, arr.ind = TRUE)
    mati <- mat_cent[kind, ]
    learnout <- learn(mat_cent[-kind, ], wavespecs, breaks)
    Extract <- learnout[[2]]
    Transform <- learnout[[3]]
    keep <- learnout[[4]]

    for (j in 1:r) {
      print(paste("= Latent Dim. =", j))
      proj_v[kind, ] <- Transform(Extract(mati, j, keep)[[1]])
      # proj_v[kind, ] <- Transform(Extract(mati, j), j)#as.matrix(t(mati))%*%as.matrix(phi_v[ ,1:j])%*%t(as.matrix(phi_v[,1:j])) # (1xp)*(pxj)*(jxp) = 1xp
      proji <- as.matrix(if (NROW(proj_v[kind, ]) > NCOL(proj_v[kind, ]) && (n < p || NROW(proj_v[kind, ]) == p)) t(proj_v[kind, ]) else proj_v[kind, ])
      matir <- if (NROW(mati) < NCOL(mati) || (n > p & NCOL(mati) == p)) as.matrix(t(mati)) else as.matrix(mati)
      rho_v[kind, j] <- sapply(
        seq.int(dim(proji)[1]),
        function(i) cor(as.matrix(proji)[i, ], as.matrix(t(matir))[i, ])^2
      )
      corM_v[j] <- cor(c(proj_v), c(mat_cent))^2
      PRESS[j] <- norm(as.matrix(proj_v - mat_cent), "F")^2
    }
  }

  # sort rows of correlation matrix
  for (s in 1:ncol(rho_v)) {
    Qrho_v[, s] <- sort(rho_v[, s])
  }

  # W Statistic?
  vline <- 0

  #* OUTPUT
  out <- list(corM_t = corM_t, rho_v = rho_v, Qrho_v = Qrho_v, vline = vline, r = r, n = n, p = p, Extract = Extract, Transform = Transform, indices = indices)
}

# see quantiles for mnist
#  can user zoom in on x axis, slider bar

# geldwt <- flf_basissel_dwt(mat=gels.data, kf=5, learn=learn_dwt
#                            , ws = list(wavelet = "db4", nlevels = 6, compress
#                                           = 0.95, ndim = 1, t = NULL, rectangular = 0, boundary = "sym")
#                            , center = TRUE, scale = FALSE
#                             , breaks=c(seq(0.8, 0.975, by=0.025), seq(0.975, 1, by=0.001)))
# # mnistdwt=gldwt
#
#
# set.seed(123)
#
# data("glaucoma_images", package = "complex.fundata"); dim(glaucoma.images)
# glauc.data <- t(matrix(glaucoma.images, 120*120, 306))
# glauc.spl <- glauc.data[sample(1:306, 18, replace=F), ]
#
# glaucdwt <- flf_basissel_dwt(mat=glauc.spl, kf=2, learn=learn_dwt
#                            , ws = list(wavelet = "db4", nlevels = c(5,5)
#                                        , compress = 1, ndim = 2, t = c(120, 120)
#                                        , rectangular = 1, boundary = c("per", "sym"))
#                            , center = TRUE, scale = FALSE
#                           , breaks=c(seq(0.8, 0.975, by=0.025), seq(0.975, 1, by=0.001)))
#
