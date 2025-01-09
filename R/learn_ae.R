#' Learning function for autoencoder (AE)
#'
#' @param Y An n times p data matrix.
#' @param k The latent feature dimension. Also known as the "bottleneck dimension" in machine learning terminology.
#' @param ae_args A list containing the following named elements to define the architecture and training of the AE: `layer_1_dim`, `link_fun`, `epochs`, `loss` and `batch_size`.
#' @return
#' @export
#'
#' @examples
#' @importFrom magrittr "%>%"

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
