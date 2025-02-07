% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/learn_dwt.2d.R
\name{learn_dwt.2d}
\alias{learn_dwt.2d}
\title{Learning Function for 2-D Discrete Wavelet Transform (DWT)}
\usage{
learn_dwt.2d(Y, p1, p2)
}
\arguments{
\item{Y}{A numeric matrix with `n` rows (observations) and `p1 * p2` columns (features),
where each row represents a flattened 2-D signal of dimensions `(p1, p2)`.}

\item{p1}{An integer specifying the number of rows in the original 2-D signal.}

\item{p2}{An integer specifying the number of columns in the original 2-D signal.}
}
\value{
A list containing:
  \itemize{
    \item `Encode`: A function to encode data into a thresholded wavelet space.
                    Takes arguments `Y` (data matrix) and `k` (number of coefficients to retain).
    \item `Decode`: A function to reconstruct data from the thresholded wavelet coefficients.
  }
}
\description{
This function applies a 2-D discrete wavelet transform (DWT) to encode and decode
image-like data (e.g., matrices or arrays) using energy-based thresholding of wavelet coefficients.
}
\details{
The function works by first reshaping the flattened input data into its original 2-D form.
Each signal is padded to a dyadic size (if necessary) to ensure compatibility with the wavelet transform.
Energy-based thresholding is then applied to retain the `k` most significant coefficients for encoding.
Reconstruction is performed using the inverse 2-D DWT.
}
