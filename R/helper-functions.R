#' Compute 1 - Squared Correlation
#'
#' Calculates 1 minus the squared correlation between observed and predicted values.
#' If `predicted` values are constant, the squared correlation is set to 0.
#'
#' @param observed Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Numeric value representing 1 - squared correlation.
#' @examples
#' get_one_minus_squared_correlation(c(1, 2, 3), c(1, 2, 2))
get_one_minus_squared_correlation <- function(observed, predicted) {
  # Handle case where predicted values are constant
  if (length(unique(predicted)) == 1) {
    warning("Predicted values constant: Setting squared correlation to 0")
    return(1)
  } else {
    return(1 - cor(observed, predicted)^2)
  }
}

#' Repeat a Column Vector
#'
#' Creates a matrix by repeating a column vector `x` `n` times.
#'
#' @param x Numeric vector to repeat.
#' @param n Integer specifying the number of repetitions.
#' @return A matrix where `x` is repeated `n` times as columns.
#' @examples
#' rep.col(c(1, 2, 3), 3)
rep.col <- function(x, n) {
  # Repeat x into a matrix with n columns
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

#' Plot Training vs Validation Loss Ratio
#'
#' Visualizes the ratio of average training loss to average validation loss as a function
#' of the number of latent features.
#'
#' @param GLaRe_output Output object from the `GLaRe()` function.
#' @return A plot showing the training-to-validation loss ratio.
#' @export
plot_train_validation_ratio <- function(GLaRe_output) {
  breaks <- GLaRe_output[["breaks"]]
  corM_t <- GLaRe_output[["corM_t"]]
  corM_v <- GLaRe_output[["corM_v"]]
  # Plot the ratio of training loss to validation loss
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

#' Plot Distribution of Validation Losses
#'
#' Creates a jittered strip chart to visualize the distribution of individual validation losses
#' for different numbers of latent features.
#'
#' @param GLaRe_output Output object from the `GLaRe()` function.
#' @param jitter Amount of jitter to add to the points. Defaults to 0.2.
#' @param cex Point size for the plot. Defaults to 1.
#' @return A strip chart of validation losses.
#' @export
distribution_plot <- function(GLaRe_output, jitter = 0.2, cex = 1) {
  Qrho_v <- GLaRe_output[["Qrho_v"]]
  breaks <- GLaRe_output[["breaks"]]
  # Plot validation losses as a strip chart
  stripchart(
    x = c(Qrho_v) ~ rep(breaks, each = nrow(Qrho_v)),
    jitter = jitter,
    method = "jitter",
    col = scales::alpha(seq_along(breaks), 0.25),
    pch = 20,
    cex = cex,
    xlab = "No. of Latent Features",
    ylab = "Individual Validation Losses",
    vertical = TRUE
  )
}

#' Plot Glaucoma Eye Data
#'
#' Visualizes glaucoma eye data as a 2D radial projection.
#'
#' @param y A 14400-dimensional numeric vector of measurements.
#' @return A ggplot2 object visualizing the eye data.
#' @export
plot_eye <- function(y) {
  # Validate input
  if (!(is.numeric(y) & length(y) == 120 * 120)) stop("y must be a numeric vector of the correct length")

  # Prepare data for plotting
  plot_df <- data.frame(
    y = y,
    theta = rep(seq(0, 360, length.out = 120), each = 120),
    phi = rep(seq(9, 24, length.out = 120), times = 120)
  )

  # Create radial plot
  ggplot2::ggplot(data = plot_df) +
    ggplot2::aes(x = theta, y = phi) +
    ggplot2::coord_radial(inner.radius = 9 / 24, expand = FALSE) +
    ggplot2::geom_tile(mapping = aes(fill = y)) +
    ggplot2::scale_fill_gradientn(colours = rainbow(8), limits = c(-0.352333, 2.020408), oob = scales::squish) +
    ggplot2::labs(fill = expression(X[i](theta, phi)), x = expression(theta), y = expression(phi))
}

#' Plot a Single MNIST Image
#'
#' Displays an individual MNIST image.
#'
#' @param y A 28x28 numeric matrix representing an MNIST image.
#' @param main Title for the plot. Defaults to NULL.
#' @return A visual plot of the MNIST image.
#' @export
plot_mnist <- function(y, main = NULL) {
  # Validate input
  stopifnot(is.matrix(y) & dim(y) == c(28, 28))
  y <- t(apply(y, 2, rev))
  # Render the image
  image(
    x = 1:28,
    y = 1:28,
    z = y,
    col = gray((0:255) / 255),
    main = main,
    ann = FALSE,
    xaxt = "n",
    yaxt = "n"
  )
}


#' Plot Eye Reconstruction
#'
#' Visualizes the original and reconstructed glaucoma eye data side by side.
#'
#' @param GLaRe_output The object returned from a call to the `GLaRe()` function.
#' @param y A 14400-dimensional numeric vector of measurements.
#' @return A ggplot2 object displaying the original and reconstructed data.
#' @export
plot_eye_reconstruction <- function(GLaRe_output, y) {
  # Extract encoding and decoding functions
  Encode <- function(Y) {
    GLaRe_output$Encode(Y)
  }
  Decode <- GLaRe_output$Decode

  # Prepare data for reconstruction
  Y <- matrix(y, nrow = 1, ncol = GLaRe_output$p)
  recon <- Decode(Ystar = Encode(Y = Y))
  recon_vec <- c(recon[1, ])

  # Generate plots for original and reconstructed data
  p1 <- plot_eye(y = y) + labs(title = "Data") + theme(plot.title = element_text(hjust = 0.5))
  p2 <- plot_eye(y = recon_vec) + labs(title = "Reconstruction") + theme(plot.title = element_text(hjust = 0.5))

  # Arrange plots side by side
  ggpubr::ggarrange(p1, p2, common.legend = TRUE)
}

#' Plot MNIST Reconstruction
#'
#' Displays the original and reconstructed MNIST image side by side.
#'
#' @param GLaRe_output The object returned from a call to the `GLaRe()` function.
#' @param y A 28x28 matrix representing an MNIST image.
#' @return A visual comparison of the original and reconstructed MNIST image.
#' @export
plot_mnist_reconstruction <- function(GLaRe_output, y) {
  # Extract encoding and decoding functions
  Encode <- GLaRe_output$Encode
  Decode <- GLaRe_output$Decode

  # Prepare data for reconstruction
  Y <- matrix(y, nrow = 1, ncol = GLaRe_output$p)
  recon <- Decode(Ystar = Encode(Y = Y))
  recon_reshape <- matrix(recon, 28, 28)
  y_mat <- matrix(y, 28, 28)

  # Plot original and reconstructed images
  par(mfrow = c(1, 2))
  plot_mnist(y = y_mat, main = "")
  title(main = "Data")
  plot_mnist(recon_reshape, main = "")
  title(main = "Reconstruction")
}

#' Plot 1D Reconstruction
#'
#' Visualizes original and reconstructed 1D data (e.g., time series).
#'
#' @param GLaRe_output The object returned from a call to the `GLaRe()` function.
#' @param Y An n x p matrix containing n observations to reconstruct.
#' @return A plot comparing the original and reconstructed 1D data.
#' @export
plot_1D_reconstruction <- function(GLaRe_output, Y) {
  # Extract encoding and decoding functions
  Encode <- GLaRe_output$Encode
  Decode <- GLaRe_output$Decode

  # Perform reconstruction
  recon <- Decode(Ystar = Encode(Y = Y))

  # Plot original and reconstructed data
  matplot(y = t(Y), type = "l", lty = 1) # Original data
  matlines(y = t(recon), type = "l", lty = 3) # Reconstructed data
  legend("topright", legend = c("Data", "Reconstruction"), lty = c(1, 3), col = 1)
}

#' Plot Gel Data
#'
#' Visualizes proteomic gels data.
#'
#' @param y A 646 x 861 matrix representing gel data.
#' @return A heatmap of the gel data.
#' @export
plot_gel <- function(y) {
  # Generate heatmap for the gel data
  image(
    x = seq(0, 1, length.out = 861),
    y = seq(0, 1, length.out = 646),
    z = y,
    col = topo.colors(10),
    xlab = NA, ylab = NA
  )
}

#' Plot Gel Reconstruction
#'
#' Visualizes the original and reconstructed proteomic gels data side by side.
#'
#' @param GLaRe_output The object returned from a call to the `GLaRe()` function.
#' @param y A 646 x 861 matrix representing gel data.
#' @return A visual comparison of the original and reconstructed gel data.
#' @export
plot_gel_reconstruction <- function(GLaRe_output, y) {
  # Prepare data for reconstruction
  Y <- matrix(y, nrow = 1, ncol = GLaRe_output$p)

  # Extract encoding and decoding functions
  Encode <- GLaRe_output$Encode
  Decode <- GLaRe_output$Decode

  # Perform reconstruction
  recon <- Decode(Ystar = Encode(Y = Y))

  # Reshape data for plotting
  Y <- matrix(y, nrow = 861, ncol = 646)
  Yhat <- matrix(recon, nrow = 861, ncol = 646)

  # Plot original and reconstructed gel data
  par(mfrow = c(1, 2))
  plot_gel(y = Y)
  title("Data")
  plot_gel(y = Yhat)
  title("Reconstruction")
}



#' Simulate Structured Data from a PCA Model
#'
#' This function generates structured data by simulating latent factors and loadings
#' from a PCA model, optionally adding Gaussian noise. The resulting data can be used
#' to test dimensionality reduction and latent representation methods.
#'
#' @param n An integer specifying the number of observations (rows) to generate.
#' @param p An integer specifying the number of features (columns) in the dataset.
#' @param k An integer specifying the number of true latent components (rank of the PCA model).
#' @param noise_sd A numeric value specifying the standard deviation of Gaussian noise
#'                 to add to the data. Defaults to 0.1.
#' @param seed An optional integer for setting the random seed to ensure reproducibility.
#'             Defaults to NULL.
#' @return A numeric matrix of dimensions `n` x `p` containing the simulated data.
#' @examples
#' # Simulate data with 100 observations, 10 features, and 3 latent components
#' data <- simulate_pca_data(n = 100, p = 10, k = 3, noise_sd = 0.1, seed = 42)
#'
#' # Plot the first two features
#' plot(data[, 1], data[, 2], main = "Scatterplot of First Two Features")
#'
#' @export
simulate_pca_data <- function(n, p, k, noise_sd = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate latent scores and loadings
  scores <- matrix(rnorm(n * k), nrow = n, ncol = k)
  loadings <- matrix(rnorm(p * k), nrow = p, ncol = k)

  # Construct data and add Gaussian noise
  data <- scores %*% t(loadings)
  data <- data + matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)

  return(data)
}
