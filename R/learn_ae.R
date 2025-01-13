#' Train an Autoencoder for Dimensionality Reduction
#'
#' This function trains an autoencoder to encode data into a lower-dimensional latent space
#' and reconstruct it back to the original space.
#'
#' @param Y A numeric matrix with `n` rows (observations) and `p` columns (variables).
#' @param k An integer specifying the latent feature dimension (bottleneck dimension).
#' @param ae_args A list containing the following hyperparameters:
#'   \itemize{
#'     \item `layer_1_dim`: Number of units in the hidden layer (default: 600).
#'     \item `link_fun`: Activation function for the output layer, either `"linear"` or `"sigmoid"` (default: `"sigmoid"`).
#'     \item `epochs`: Number of training epochs (default: 100).
#'     \item `loss`: Loss function, either `"mean_squared_error"` or `"binary_crossentropy"` (default: `"mean_squared_error"`).
#'     \item `batch_size`: Mini-batch size for training (default: 16).
#'   }
#' @return A list containing:
#'   \itemize{
#'     \item `Encode`: A function to encode data into the latent space.
#'     \item `Decode`: A function to reconstruct data from the latent space.
#'   }
#' @examples
#' # Example: Simulated Data from a PCA Model
#' library(GLarE)
#' set.seed(1996)
#'
#' # Simulate data from a PCA model
#' n <- 100 # Number of observations
#' p <- 10 # Number of variables
#' k <- 3 # True number of latent components
#' loadings <- matrix(rnorm(p * k), nrow = p, ncol = k) # Loadings matrix
#' scores <- matrix(rnorm(n * k), nrow = n, ncol = k) # Factor scores
#' Y <- scores %*% t(loadings) + matrix(rnorm(n * p, mean = 0, sd = 0.1), nrow = n, ncol = p)
#' # Train autoencoder
#' ae_model <- learn_ae(Y, k = 5, ae_args = list(epochs = 50, layer_1_dim = 100, link_fun = "linear"))
#' encoded_data <- ae_model$Encode(Y)
#' reconstructed_data <- ae_model$Decode(encoded_data)
#' @importFrom magrittr "%>%"
#' @export

learn_ae <- function(Y, k, ae_args) {
  # unpack arguments:
  layer_1_dim <- ifelse(is.null(ae_args[["layer_1_dim"]]), 600, ae_args[["layer_1_dim"]])
  link_fun <- ifelse(is.null(ae_args[["link_fun"]]), "sigmoid", ae_args[["link_fun"]])
  epochs <- ifelse(is.null(ae_args[["epochs"]]), 100, ae_args[["epochs"]])
  loss <- ifelse(is.null(ae_args[["loss"]]), "mean_squared_error", ae_args[["loss"]])
  batch_size <- ifelse(is.null(ae_args[["batch_size"]]), 16, ae_args[["batch_size"]])


  if (!(link_fun %in% c("sigmoid", "linear"))) stop("Link function must be either linear or sigmoid.")
  if (!(loss %in% c("mean_squared_error", "binary_crossentropy"))) stop("Loss must be mean_squared_error or binary_crossentropy")

  p <- ncol(Y)

  # Define the encoder
  encoder <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = layer_1_dim, activation = "relu", input_shape = p) %>%
    keras::layer_dense(units = k, activation = "linear", name = "bottleneck")

  # Define the decoder
  decoder <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = layer_1_dim, activation = "relu", input_shape = k) %>%
    keras::layer_dense(units = p, activation = link_fun)

  # Connect them to create the autoencoder
  autoencoder <- keras::keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder %>% keras::compile(optimizer = "adam", loss = loss)
  autoencoder %>% keras::fit(Y,
    Y,
    epochs = epochs,
    batch_size = batch_size,
    verbose = FALSE
  )

  Encode <- function(Y) {
    predict(encoder, Y, verbose = FALSE)
  }

  Decode <- function(Ystar) {
    predict(decoder, Ystar, verbose = FALSE)
  }

  list(Encode = Encode, Decode = Decode)
}
