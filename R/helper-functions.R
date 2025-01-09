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
#' @param GLaRe_output the object returned from a call to the `GLaRE()` function.
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 "ggplot"
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
#' @param GLaRe_output the object returned from a call to the `GLaRE()` function.
#'
#' @return
#' @export
#'
#' @examples
distribution_plot <- function(GLaRe_output, jitter = 0.2, cex = 1) {
  Qrho_v <- GLaRe_output[["Qrho_v"]]
  breaks <- GLaRe_output[["breaks"]]
  stripchart(
    x = c(Qrho_v) ~ rep(breaks, each = nrow(Qrho_v)),
    jitter = 0.2,
    method = "jitter",
    col = scales::alpha(seq_along(breaks), 0.25),
    pch = 20,
    cex = cex,
    xlab = "No. of Latent Features",
    ylab = "Individual Validation Losses",
    vertical = TRUE
  )
}


#' Plot the eye glaucoma data.
#'
#' @param y a 14400-dimensional vector comprising measurements of MPS along `theta = rep(seq(0, 360, length.out = 120)`, each = 120) and `phi = rep(seq(9, 24, length.out = 120), times = 120)`.
#'
#' @return
#' @export
#'
#' @examples
#' @import ggplot2
plot_eye <- function(y) {
  if (!(is.numeric(y) & length(y) == 120 * 120)) stop("y must be a numeric vector of measurements on theta = rep(seq(0, 360, length.out = 120), each = 120) and phi = rep(seq(9, 24, length.out = 120), times = 120)")
  plot_df <- data.frame(
    y = y,
    theta = rep(seq(0, 360, length.out = 120), each = 120),
    phi = rep(seq(9, 24, length.out = 120), times = 120)
  )

  p <- ggplot2::ggplot(data = plot_df) +
    ggplot2::aes(x = theta, y = phi) +
    ggplot2::coord_radial(inner.radius = 9 / 24, expand = FALSE) +
    ggplot2::geom_tile(mapping = aes(fill = y)) +
    ggplot2::scale_fill_gradientn(colours = rainbow(8), limits = c(-0.352333, 2.020408), oob = scales::squish) +
    ggplot2::labs(fill = expression(X[i](theta, phi)), x = expression(theta), y = expression(phi))
  p
}


#' Helper function to plot a single image from the mnist data.
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
plot_mnist <- function(y, main = NULL) {
  stopifnot(is.matrix(y) & dim(y) == c(28, 28))
  y <- t(apply(y, c(2), rev))
  image(
    x = 1:28,
    y = 1:28,
    z = y,
    col = gray((0:255) / 255),
    main = main,
    ann = F,
    xaxt = "n",
    yaxt = "n"
  )
}


#' Plot a reconstruction of the glaucoma data.
#'
#' @param GLaRe_output the object returned from a call to the `GLaRE()` function.
#' @param y a 14400-dimensional vector comprising measurements of MPS along `theta = rep(seq(0, 360, length.out = 120)`, each = 120) and `phi = rep(seq(9, 24, length.out = 120), times = 120)`.
#'
#' @return
#' @export
#'
#' @examples
plot_eye_reconstruction <- function(GLaRe_output, y) {
  Encode <- function(Y) {
    GLaRe_output$Encode(Y)
  }
  Decode <- GLaRe_output$Decode
  Y <- matrix(y, nrow = 1, ncol = GLaRe_output$p)
  recon <- Decode(Ystar = Encode(Y = Y))
  recon_vec <- c(recon[1, ])
  p1 <- plot_eye(y = y) + labs(title = "Data") + theme(plot.title = element_text(hjust = 0.5))
  p2 <- plot_eye(y = recon_vec) + labs(title = "Reconstruction") + theme(plot.title = element_text(hjust = 0.5))
  ggpubr::ggarrange(p1, p2, common.legend = TRUE)
}


#' Plot a reconstruction of the mnist data.
#'
#' @param GLaRe_output the object returned from a call to the `GLaRE()` function.
#' @param y a 28 times 28 matrix representing an image from the MNIST dataset.
#'
#' @return
#' @export
#'
#' @examples
plot_mnist_reconstruction <- function(GLaRe_output, y) {
  Encode <- GLaRe_output$Encode
  Decode <- GLaRe_output$Decode
  Y <- matrix(y, nrow = 1, ncol = GLaRe_output$p)
  recon <- Decode(Ystar = Encode(Y = Y))
  recon_reshape <- matrix(recon, 28, 28)
  y_mat <- matrix(y, 28, 28)
  par(mfrow = c(1, 2))
  plot_mnist(y = y_mat, main = "")
  title(main = "Data")
  plot_mnist(recon_reshape, main = "")
  title(main = "Reconstruction")
}

#' Plot a reconstruction of 1-D (e.g., time series) data.
#'
#' @param GLaRe_output the object returned from a call to the `GLaRE()` function.
#' @param Y An n times p matrix containing the n observations for which you want their reconstruction displayed.
#' @return
#' @export
#'
#' @examples
plot_1D_reconstruction <- function(GLaRe_output, Y) {
  Encode <- GLaRe_output$Encode
  Decode <- GLaRe_output$Decode
  recon <- Decode(Ystar = Encode(Y = Y))
  matplot(y = t(Y), type = "l", lty = 1)
  matlines(y = t(recon), type = "l", lty = 3)
  legend("topright", legend = c("Data", "Reconstruction"), lty = c(1, 3), col = 1)
}

plot_gel <- function(y) {
  image(x = seq(0, 1, length.out = 861),
        y = seq(0, 1, length.out = 646),
        z = y,
        col = topo.colors(10),
        xlab = NA, ylab = NA)
}

#' Plot a reconstruction of gels data.
#'
#' @param GLaRe_output the object returned from a call to the `GLaRE()` function.
#' @param Y An n times p matrix containing the n observations for which you want their reconstruction displayed.
#' @return
#' @export
#'
#' @examples
plot_gel_reconstruction <- function(GLaRe_output, y) {

  Y <- matrix(y, nrow = 1, ncol = GLaRe_output$p)


  Encode <- GLaRe_output$Encode
  Decode <- GLaRe_output$Decode

  recon <- Decode(Ystar = Encode(Y = Y))

  Y <- matrix(Y, nrow = 861, ncol = 646)
  Yhat <- matrix(recon, nrow = 861, ncol = 646)

  par(mfrow = c(1, 2))
  plot_gel(y = Y)
  title("Data")
  plot_gel(y = Yhat)
  title("Reconstruction")
}
