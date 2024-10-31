get_one_minus_squared_correlation <- function(observed, predicted) {
  # check predicted
  if (length(unique(predicted)) == 1) {
    warning("Predicted values constant: Setting squared correlation to 0")
    return(1)
  } else {
    return(1 - cor(observed, predicted)^2)
  }
}

rep.col <- function(x, n) { # this function repeats a column x n times
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

#' Plot the ratio of the average training loss to the average validation loss.
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
plot_train_validation_ratio <- function(GLaRe_output) {
  breaks <- GLaRe_output[["breaks"]]
  corM_t <- GLaRe_output[["corM_t"]]
  corM_v <- GLaRe_output[["corM_v"]]
  plot(
    x = breaks,
    y = corM_t / corM_v,
    type = "b",
    ylim = c(0, 1),
    pch = 20,
    xlab = "No. of Latent Features",
    ylab = "Ratio of Training to Validation Loss"
  )
}

#' Plot the full distribution of individual validation losses.
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
distribution_plot <- function(GLaRe_output, jitter = 0.2, cex = 1) {
  Qrho_v <- GLaRe_output[["Qrho_v"]]
  breaks <- GLaRe_output[["breaks"]]
  stripchart(x = c(Qrho_v) ~ rep(breaks, each = nrow(Qrho_v)),
             jitter = 0.2,
             method="jitter",
             col = scales::alpha(seq_along(breaks), 0.25),
             pch = 20,
             cex = cex,
             xlab = "No. of Latent Features",
             ylab = "Individual Validation Losses",
             vertical = TRUE)
}
