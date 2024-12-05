#' Pad a signal to dyadic length for 1-D DWT with 0's each side.
#'
#' @param Y an N times T data matrix containing N replicate observations of a T-dimensional signal in its rows.
prepare_pad_dwt <- function(Y) {
  p <- ncol(Y)
  log2ppad <- ceiling(log2(p))
  ppad <- 2^log2ppad
  ppad_extra <- ppad - p
  ppad_left <- ceiling(ppad_extra / 2)
  ppad_right <- floor(ppad_extra / 2)
  if (p != ppad) {
    Ypad <- cbind(matrix(0, nrow(Y), ncol = ppad_left), Y, matrix(0, nrow(Y), ncol = ppad_right))
  } else if (p == ppad) {
    Ypad <- Y
  }
  list(Ypad = Ypad, ppad = ppad, log2ppad = log2ppad, ppad_left = ppad_left, ppad_right = ppad_right)
}

#' Apply DWT to a 1-D signal and return the wavelet coefficients as a list.
#' @param y a vector conatining the 1-D signal, must be of dyadic length.
waveslim_dwt <- function(y) {
  dwt_obj <- waveslim::dwt(x = y, n.levels = log2(length(y)), wf = "la8", boundary = "periodic")
  lapply(dwt_obj, function(x) {
    x
  })
}

#' Apply DWT to a 1-D signal and return the wavelet coefficients as a vector.
#' @param y a vector conatining the 1-D signal, must be of dyadic length.
waveslim_dwt_to_vec <- function(y) {
  dwt_list <- waveslim_dwt(y)
  unlist(dwt_list)
}

#' Apply DWT to a a matrix of n 1-D signals and return the wavelet coefficients as a matrix.
#' @param Y an N times T data matrix containing N replicate observations of a T-dimensional signal in its rows. T must be dyadic length.
waveslim_dwt_to_mat <- function(Y) {
  t(apply(Y, 1, waveslim_dwt_to_vec))
}

#' Compute relative energy for a general wavelet coefficient matrix.
#'
#' @param D An N x T matrix of wavelet coefficients.
#' @return En A N x T matrix containing the relative energy, such that row i and
#' column k entry contains the cumulative proportion of energy explained by this
#' coefficient and all coefficients greater in absolute value than it for subject i.
Energy <- function(D) {
  En <- D # energy matrix has same dimensions of coef matrix.
  total_Energy <- vector("numeric", length = nrow(D)) # each subject has a total energy
  for (i in 1:nrow(D)) {
    # indexing subject (row) i....
    ix <- order(abs(D[i, ]), decreasing = TRUE) # indices of wavelet coefficients, ordered from largest to smallest abs value
    Csort <- D[i, ix] # vector of wavelet coefficients ordered according size (largest to smallest)
    total_Energy[i] <- sum(Csort^2) # total energy = sum of all squared coefs (could have computed before order)
    energy <- cumsum(as.numeric(Csort^2)) / total_Energy[i] # ordered relative energy vector: kth entry contains cumulative prop of energy explained by k largest wavelet coefs.
    En[i, ix] <- energy # vector such that entry k = cumulative proportion of energy explained by this coefficient and all coefficients greater in absolute value than it.
  }
  En
}

#' Compute scree from relative energy matrix.
#'
#' @param D An N x T matrix of wavelet coefficients.
#' @return En A N x T matrix containing the relative energy, such that row i and
#' column k entry contains the cumulative proportion of energy explained by this
#' coefficient and all coefficients greater in absolute value than it.
get_Energy_scree <- function(D) {
  Dcomp <- Energy(D) # get relative energy matrix.
  scree <- apply(Dcomp, 2, mean) # calculate mean energy for each coefficient.
}

threshold_fun <- function(D, scree, k) {
  screeix <- order(scree, decreasing = FALSE) # order coefs by mean energy explained (smallest to largest)
  # intuition here is that if rel energy = 1, then all energy has been explained before this coef.
  Dthresh_inds <- screeix[-c(1:k)] # remove the k largest scree inds
  D0 <- D
  D0[, Dthresh_inds] <- 0 # set the remainder (i.e., not the k largest scree inds) to zero.
  D0 # return thresholded matrix.
}

idwt_vec <- function(d, ppad, ppad_left, ppad_right) {
  dwt_obj_refill <- waveslim::dwt(rep(0, ppad), n.levels = log2(ppad), wf = "la8", boundary = "periodic")
  jinds <- seq_along(dwt_obj_refill[[1]])
  dwt_obj_refill[[1]] <- d[jinds]
  for (j in 2:length(dwt_obj_refill)) {
    offset_j <- sum(sapply(dwt_obj_refill[1:(j - 1)], length))
    jinds <- offset_j + seq_along(dwt_obj_refill[[j]])
    dwt_obj_refill[[j]] <- d[jinds]
  }
  waveslim::idwt(y = dwt_obj_refill)[(ppad_left + 1):(ppad - ppad_right)]
}

idwt_mat <- function(D, ppad, ppad_left, ppad_right) {
  t(apply(D, 1, idwt_vec, ppad = ppad, ppad_left = ppad_left, ppad_right = ppad_right))
}

#' Learning function for 1-D DWT.
#'
#' @param Y an n times p data matrix.
#'
learn_dwt <- function(Y) {
  n <- nrow(Y)
  p <- ncol(Y)
  Ypad_obj <- prepare_pad_dwt(Y = Y)
  Ypad <- Ypad_obj$Ypad
  ppad <- Ypad_obj$ppad
  ppad_left <- Ypad_obj$ppad_left
  ppad_right <- Ypad_obj$ppad_right
  D_train <- waveslim_dwt_to_mat(Y = Ypad)
  scree <- get_Energy_scree(D = D_train)

  Encode <- function(Y, k) {
    Ypad <- prepare_pad_dwt(Y = Y)$Ypad
    D <- waveslim_dwt_to_mat(Y = Ypad)
    threshold_fun(D = D, scree = scree, k = k)
  }

  Decode <- function(Ystar) {
    idwt_mat(D = Ystar, ppad = ppad, ppad_left = ppad_left, ppad_right = ppad_right)
  }

  list(Encode = Encode, Decode = Decode)
}
