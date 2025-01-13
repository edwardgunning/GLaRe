transform_correlation_output <- function(out_basissel, cvqlines, cutoff_criterion) {
  if (!(cvqlines >= 0 & cvqlines <= 1)) stop("cvqlines must be in [0, 1]")
  if (!(is.list(out_basissel))) stop("out_basissel must be the named list output")
  if (!all(c("corM_t", "rho_v") %in% names(out_basissel))) stop("out_basissel object must contain elements corM_t and rho_v")
  cor_df <- data.frame(
    meansqcor_t = out_basissel[["corM_t"]],
    minsqcor_cv = apply(out_basissel[["rho_v"]], 2, min),
    meansqcor_cv = out_basissel[["corM_v"]],
    maxsqcor_cv = apply(out_basissel[["rho_v"]], 2, max),
    medsqcor_cv = apply(out_basissel[["rho_v"]], 2, function(x) quantile(x, 0.5, na.rm = TRUE)),
    qchoice_cv = apply(out_basissel[["rho_v"]], 2, function(x) quantile(x, cvqlines, na.rm = TRUE)),
    cutoff_criterion_cv = apply(out_basissel[["rho_v"]], 2, function(x) quantile(x, cutoff_criterion, na.rm = TRUE))
  )
  cor_df
}

summary_correlation_plot <- function(out_basisel, cvqlines, cutoff_criterion, r, q, breaks, method_name, qd, tolerance_level) {
  correlation_df <- transform_correlation_output(out_basisel, cvqlines, cutoff_criterion)
  plot(
    x = breaks,
    y = correlation_df[, "meansqcor_t"],
    type = "b",
    pch = 20,
    col = "green",
    lwd = 3,
    lty = 1,
    panel.first = c(abline(h = 0, lty = 1, col = "black"), abline(h = 1, lty = 1, col = "black")),
    xlab = "No. of Latent Features",
    ylab = expression(paste("Loss: 1 - Squared Correlation (", 1 - R^2, ")", sep = "")),
    xlim = range(breaks),
    ylim = c(0, 1),
    main = paste("Latent Feature Representation \n Summary:", method_name)
  )
  lines(x = breaks[seq_len(r)], correlation_df[seq_len(r), "minsqcor_cv"], col = "royalblue", lwd = 2, type = "b", pch = 20)
  lines(x = breaks[seq_len(r)], correlation_df[seq_len(r), "meansqcor_cv"], col = "goldenrod", lwd = 2, type = "b", pch = 20)
  lines(x = breaks[seq_len(r)], correlation_df[seq_len(r), "maxsqcor_cv"], col = "red3", lwd = 2, type = "b", pch = 20)
  if (cvqlines == 0.5) {
    lines(x = breaks[seq_len(r)], correlation_df[seq_len(r), "medsqcor_cv"], col = "purple", lwd = 2, type = "b", pch = 20)
  } else if (cvqlines != 0.5) {
    lines(x = breaks[seq_len(r)], correlation_df[seq_len(r), "qchoice_cv"], col = "purple", lwd = 2, type = "b", pch = 20)
  }
  lines(x = breaks[seq_len(r)], correlation_df[seq_len(r), "cutoff_criterion_cv"], col = "grey", lwd = 2, lty = 2, type = "b", pch = 20)

  if (!is.na(qd)) {
    abline(v = qd, lty = 2, col = "grey")
    abline(h = tolerance_level, lty = 2, col = "grey")
    axis(side = 1, at = c(qd), labels = paste0("qd = ", qd), col = "darkgrey", font = 4, lwd = 3, padj = 1.2)
    axis(side = 2, at = c(tolerance_level), labels = paste0("ε = ", tolerance_level), col = "darkgrey", font = 4, lwd = 3, padj = 1.2)
  }


  legend("topright",
    legend = c(
      "CV Min Loss",
      "CV Mean Loss",
      paste("CV Percentile =", cvqlines, "Loss"),
      "CV Max Loss",
      "Training Mean Loss",
      paste("Cut-Off Criterion = ", cutoff_criterion, "Loss")
    ),
    col = c("blue", "goldenrod", "purple", "red3", "green", "grey"),
    lty = c(1, 1, 1, 1, 1, 2),
    lwd = c(2, 2, 2, 2, 2, 2, 2),
    bg = "white"
  )
}


#' Assess Losslessness of a Latent Feature Representation Method using k-Fold Cross-Validation
#'
#' The `GLaRe` function evaluates the quality of latent feature representation methods by estimating
#' the distribution of information loss through k-fold cross-validation. It identifies the qualifying
#' dimension (qd) where a user-defined quantile of the loss distribution meets a tolerance level,
#' providing a compact and interpretable representation. Visualizations such as scree plots and
#' heatmaps are generated for interpretability.
#'
#' @param mat An n-by-p data matrix where `n` represents the number of observations, and `p` represents
#'        the number of features.
#' @param latent_dim_from An integer specifying the starting point of the latent feature dimension
#'        range. Defaults to 1.
#' @param latent_dim_to An integer specifying the ending point of the latent feature dimension
#'        range. Defaults to `min(ncol(mat) - 1, nrow(mat) - 1)`.
#' @param latent_dim_by An integer defining the step size for increments in the latent feature dimension
#'        range. Defaults to 1.
#' @param learn A string specifying the latent feature representation method. Options are:
#'        \itemize{
#'          \item `"pca"`: Principal Component Analysis.
#'          \item `"dwt"`: Discrete Wavelet Transform.
#'          \item `"dwt.2d"`: Two-dimensional Discrete Wavelet Transform.
#'          \item `"ae"`: Autoencoder.
#'          \item `"user"`: User-defined custom method (requires `learn_function` parameter).
#'        }
#'        Defaults to `"pca"`.
#' @param method_name A string for the name of the method to feature on plot titles. Defaults to
#'        `toupper(learn)`.
#' @param kf An integer specifying the number of folds for k-fold cross-validation. Defaults to 5.
#' @param cvqlines A numeric value between 0 and 1 specifying the quantile of the cross-validated
#'        loss distribution to display on the scree plot. Defaults to 0.9.
#' @param ae_args A list of parameters for autoencoder training, used only if `learn = "ae"`. Must include:
#'        \itemize{
#'          \item `layer_1_dim`: Number of nodes in the hidden layer.
#'          \item `link_fun`: Activation function for the autoencoder.
#'          \item `epochs`: Number of training epochs.
#'          \item `loss`: Loss function to optimize.
#'          \item `batch_size`: Mini-batch size for training.
#'        }
#' @param tolerance_level A numeric value specifying the maximum allowable loss for a specified quantile
#'        of observations. Defaults to 0.05.
#' @param cutoff_criterion A numeric value (between 0 and 1) defining the quantile of observations for
#'        which the tolerance level must be met. Defaults to 0.95.
#' @param learn_function A custom function for encoding and decoding, required if `learn = "user"`.
#'        The function must take arguments `Y` (data matrix) and `k` (latent dimension) and return a
#'        list containing `Encode` and `Decode` functions.
#' @param verbose Logical; if `TRUE`, progress messages are printed during training and evaluation.
#'        Defaults to `TRUE`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `qd`: The qualifying dimension where the loss meets the specified tolerance level.
#'   \item `Encode`: Encoding function for the latent representation (if `qd` is found).
#'   \item `Decode`: Decoding function for reconstructing data from the latent representation.
#'   \item `heatmap`: A Plotly heatmap of information loss across latent dimensions and quantiles.
#' }
#' @export
#' @importFrom magrittr "%>%"
#'
#' @examples
#' data(glaucoma_data)
#' result <- GLaRe(
#'   mat = as.matrix(glaucoma_data),
#'   latent_dim_from = 1,
#'   latent_dim_to = 50,
#'   learn = "pca",
#'   kf = 5,
#'   tolerance_level = 0.05,
#'   cutoff_criterion = 0.95,
#'   verbose = TRUE
#' )
#'
#' @seealso \code{\link{summary_correlation_plot}}, \code{\link{plot_train_validation_ratio}}

GLaRe <- function(
    mat = as.matrix(glaucoma_data),
    latent_dim_from = 1,
    latent_dim_to = min(ncol(mat), nrow(mat) - 1),
    latent_dim_by = 1,
    learn = "pca",
    method_name = toupper(learn),
    kf = 5,
    cvqlines = 0.9,
    ae_args = list(),
    tolerance_level = 0.05,
    cutoff_criterion = 0.95,
    learn_function = NULL,
    verbose = TRUE) {

  # start:
  if(verbose) {
    print(paste("*** Learning Method:", learn, "***"))
  }
  # apply training and cross validation for basis selection:
  out <- flf_basissel(mat = mat, learn = learn, kf = kf, latent_dim_from = latent_dim_from, latent_dim_to = latent_dim_to, latent_dim_by = latent_dim_by, ae_args = ae_args, learn_function = learn_function, verbose = verbose)
  n <- out[["n"]]
  q <- out[["q"]]
  r <- out[["r"]]
  breaks <- out[["breaks"]]

  # Qualifying Criterion and Add to plot: -----------------------------------
  cutoff_criterion_quantiles <- apply(out[["rho_v"]][, seq_len(r), drop = FALSE], 2, function(x) quantile(x, cutoff_criterion, na.rm = TRUE))
  if (!any(cutoff_criterion_quantiles <= tolerance_level)) {
    warning("No qualifying criterion found, try adjusting parameters.")
    qd <- NA
  } else {
    index <- min(which(cutoff_criterion_quantiles <= tolerance_level))
    qd <- breaks[index]
  }

  # Loss-lessness (Scree) Plot:  ---------------------------------------------
  summary_correlation_plot(out_basisel = out,
                           cvqlines = cvqlines,
                           cutoff_criterion = cutoff_criterion,
                           r = r,
                           q = q,
                           breaks = breaks,
                           method_name = method_name,
                           qd = qd,
                           tolerance_level = tolerance_level)


  # Heatmap: ----------------------------------------------------------------
 heat_map <- plotly::plot_ly(
    x = breaks,
    y = seq(0, 1, length.out = n),
    z = log10(out$Qrho_v),
    type = "heatmap",
    colorscale = list(
      list(0, "green"), # Start color (was red)
      list(0.5, "yellow"), # Middle color
      list(1, "red") # End color (was green)
    )
  ) %>%
    plotly::layout(
      # title = "Squared Correlation Heatmap",
      xaxis = list(title = "No. of Latent Features"),
      yaxis = list(title = "Quantile of Observations"),
      legend = list(orientation = "h")
    ) %>%
    plotly::colorbar(
      title = "Loss",
      tickvals = log10(10^seq(0, -10, by = -1)),
      ticktext = 10^seq(0, -10, by = -1),
      side = "bottom"
    )



  # Return Decode and Encode: -----------------------------------------------
  # Only if qualifying criterion is met:
  if(!is.na(qd)) {
    if(verbose) {
      print("Final training of Model at qualifying criterion:")
    }

    if(learn == "dwt.2d") {
    arr <- mat
    arr_dims <- dim(arr)
    n <- arr_dims[1]
    p1 <- arr_dims[2]
    p2 <- arr_dims[3]
    p <- p1 * p2
    mat <- matrix(arr, nrow = n, ncol = p)
    }

    if (learn == "pca") {
      learn_function <- learn_pca
    } else if (learn == "dwt") {
      learn_function <- learn_dwt
    } else if (learn == "dwt.2d") {
      learn_function <- function(Y) {
        learn_dwt.2d(Y = Y, p1 = p1, p2 = p2)
      }
    } else if (learn == "ae") {
      learn_function <- function(Y, k) {
        learn_ae(Y = Y, k = k, ae_args = ae_args)
      }
    }

    if (learn %in% c("pca", "dwt", "dwt.2d")) {
      learn_final <- learn_function(mat)
      out[["Encode"]] <- function(Y) {
        learn_final[["Encode"]](Y = Y, k = qd)
      }
    } else {
      learn_final <- learn_function(Y = mat, k = qd)
      out[["Encode"]] <- learn_final[["Encode"]]
    }

    out[["Decode"]] <- learn_final[["Decode"]]

  } else if(is.na(qd)){
    out[["Encode"]] <- NA
    out[["Decode"]] <- NA
  }

  out$qd <- qd
  out$heatmap <- heat_map
  # final:
  out
}
