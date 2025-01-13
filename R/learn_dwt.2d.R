prepare_pad_dwt.2d <- function(Y) {
  if (!length(dim(Y)) == 3) {
    stop("Y must be a 3-dimensional array")
  }

  Y_dims <- dim(Y)
  n <- Y_dims[1]
  p <- Y_dims[2]
  q <- Y_dims[3]

  log2ppad <- ceiling(log2(p))
  ppad <- 2^log2ppad
  ppad_extra <- ppad - p
  ppad_left <- ceiling(ppad_extra / 2)
  ppad_right <- floor(ppad_extra / 2)

  log2qpad <- ceiling(log2(q))
  qpad <- 2^log2qpad
  qpad_extra <- qpad - q
  qpad_left <- ceiling(qpad_extra / 2)
  qpad_right <- floor(qpad_extra / 2)

  if (!all(c(p, q) == c(ppad, qpad))) {
    Ypad <- array(data = 0, dim = c(n, ppad, qpad))
    Ypad[, (ppad_left + 1):(ppad - ppad_right), (qpad_left + 1):(qpad - qpad_right)] <- Y
  } else if (c(p, q) == c(ppad, qpad)) {
    Ypad <- Y
  }

  list(
    Ypad = Ypad,
    ppad = ppad, log2ppad = log2ppad, ppad_left = ppad_left, ppad_right = ppad_right,
    qpad = qpad, log2qpad = log2qpad, qpad_left = qpad_left, qpad_right = qpad_right
  )
}


waveslim_dwt.2d <- function(y) {
  p <- nrow(y)
  q <- ncol(y)
  n.levels <- log2(min(p, q))
  dwt.2d_obj <- waveslim::dwt.2d(x = y, J = n.levels, wf = "la8", boundary = "periodic")
  lapply(dwt.2d_obj, function(x) {
    x
  })
}

waveslim_dwt.2d_to_vec <- function(y) {
  dwt.2d_list <- waveslim_dwt.2d(y)
  unlist(dwt.2d_list)
}

waveslim_dwt.2d_to_mat <- function(Y) {
  t(apply(Y, 1, waveslim_dwt.2d_to_vec))
}



idwt.2d_vec <- function(d, ppad, ppad_left, ppad_right, qpad, qpad_left, qpad_right) {
  Ynull <- matrix(data = 0, nrow = ppad, ncol = qpad)
  n.levels <- log2(min(ppad, qpad))
  dwt.2d_obj_refill <- waveslim::dwt.2d(Ynull, J = n.levels, wf = "la8", boundary = "periodic")
  jinds <- seq_along(dwt.2d_obj_refill[[1]])
  nrow_j <- nrow(dwt.2d_obj_refill[[1]])
  ncol_j <- ncol(dwt.2d_obj_refill[[1]])
  dwt.2d_obj_refill[[1]] <- matrix(d[jinds], nrow = nrow_j, ncol = ncol_j)
  for (j in 2:length(dwt.2d_obj_refill)) {
    offset_j <- sum(sapply(dwt.2d_obj_refill[1:(j - 1)], length))
    nrow_j <- nrow(dwt.2d_obj_refill[[j]])
    ncol_j <- ncol(dwt.2d_obj_refill[[j]])
    jinds <- offset_j + seq_along(dwt.2d_obj_refill[[j]])
    dwt.2d_obj_refill[[j]] <- matrix(d[jinds], nrow = nrow_j, ncol = ncol_j)
  }
  waveslim::idwt.2d(y = dwt.2d_obj_refill)[(ppad_left + 1):(ppad - ppad_right), (qpad_left + 1):(qpad - qpad_right)]
}

idwt.2d_array <- function(D, ppad, ppad_left, ppad_right, qpad, qpad_left, qpad_right, p, q) {
  n <- nrow(D)
  idwt_array <- array(data = NA, dim = c(n, p, q))
  for (i in seq_len(nrow(D))) {
    idwt_array[i, , ] <- idwt.2d_vec(d = D[i, ], ppad = ppad, ppad_left = ppad_left, ppad_right = ppad_right, qpad = qpad, qpad_left = qpad_left, qpad_right = qpad_right)
  }
  idwt_array
}


#' Learning Function for 2-D Discrete Wavelet Transform (DWT)
#'
#' This function applies a 2-D discrete wavelet transform (DWT) to encode and decode
#' image-like data (e.g., matrices or arrays) using energy-based thresholding of wavelet coefficients.
#'
#' @param Y A numeric matrix with `n` rows (observations) and `p1 * p2` columns (features),
#'          where each row represents a flattened 2-D signal of dimensions `(p1, p2)`.
#' @param p1 An integer specifying the number of rows in the original 2-D signal.
#' @param p2 An integer specifying the number of columns in the original 2-D signal.
#' @return A list containing:
#'   \itemize{
#'     \item `Encode`: A function to encode data into a thresholded wavelet space.
#'                     Takes arguments `Y` (data matrix) and `k` (number of coefficients to retain).
#'     \item `Decode`: A function to reconstruct data from the thresholded wavelet coefficients.
#'   }
#' @details
#' The function works by first reshaping the flattened input data into its original 2-D form.
#' Each signal is padded to a dyadic size (if necessary) to ensure compatibility with the wavelet transform.
#' Energy-based thresholding is then applied to retain the `k` most significant coefficients for encoding.
#' Reconstruction is performed using the inverse 2-D DWT.
#'
#' @examples
#'
#' @importFrom waveslim dwt.2d idwt.2d
#' @export
learn_dwt.2d <- function(Y, p1, p2) {

  stopifnot(ncol(Y) == p1 * p2)
  n <- nrow(Y)
  Y <- array(Y, dim = c(n, p1, p2))

  Ypad_obj <- prepare_pad_dwt.2d(Y = Y)
  Ypad <- Ypad_obj$Ypad
  ppad <- Ypad_obj$ppad
  qpad <- Ypad_obj$qpad
  ppad_left <- Ypad_obj$ppad_left
  ppad_right <- Ypad_obj$ppad_right
  qpad_left <- Ypad_obj$qpad_left
  qpad_right <- Ypad_obj$qpad_right
  D_train <- waveslim_dwt.2d_to_mat(Y = Ypad)
  scree <- get_Energy_scree(D = D_train)

  Encode <- function(Y, k) {
    Y <- array(Y, dim = c(nrow(Y), p1, p2))
    Ypad <- prepare_pad_dwt.2d(Y = Y)$Ypad
    D <- waveslim_dwt.2d_to_mat(Y = Ypad)
    threshold_fun(D = D, scree = scree, k = k)
  }

  Decode <- function(Ystar) {
    arr_return <- idwt.2d_array(
      D = Ystar,
      ppad = ppad,
      ppad_left = ppad_left,
      ppad_right = ppad_right,
      qpad = qpad,
      qpad_left = qpad_left,
      qpad_right = qpad_right,
      p = p1,
      q = p2
    )
    matrix(arr_return, nrow = nrow(Ystar), ncol = p1 * p2)
  }
  list(Encode = Encode, Decode = Decode)
}
