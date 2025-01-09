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

summary_correlation_plot <- function(out_basisel, cvqlines, cutoff_criterion, r, q, breaks, method_name, qc, tolerance_level) {
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

  if (!is.na(qc)) {
    abline(v = qc, lty = 2, col = "grey")
    abline(h = tolerance_level, lty = 2, col = "grey")
    axis(side = 1, at = c(qc), labels = paste0("qc = ", qc), col = "darkgrey", font = 4, lwd = 3, padj = 1.2)
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


#' Assess losslessness of a latent feature representation method using k-fold cross-validation.
#'
#' @param mat An n-by-p data matrix (n observations, p variables).
#' @param latent_dim_from An integer used to define the start of the sequence of the number of latent features to be defined as `seq(from = latent_dim_from, to = latent_dim_to, by = latent_dim_by)`, defaults to 1.
#' @param latent_dim_to An integer used to define the end of the sequence of the number of latent features to be defined as `seq(from = latent_dim_from, to = latent_dim_to, by = latent_dim_by)` defaults to `min(ncol(mat) - 1, nrow(mat) - 1)`.
#' @param latent_dim_by An integer used to define the increment of the sequence of the number of latent features to be defined as `seq(from = latent_dim_from, to = latent_dim_to, by = latent_dim_by)`, defaults to 1.
#' @param learn The latent feature representation method chosen, one of c("pca", "dwt", "dwt.2d", "ae", "user"). Defaults to "pca" for principal component analysis (PCA).
#' @param kf An integer defining the number of folds for the k-fold cross-validation.
#' @param method_name The name of the method to be featured on the plot title. Defaults to `toupper(learn)`.
#' @param cvqlines The user-specified quantile of the cross-validated loss distribution to display on the plot.
#' @param ae_args Only to be specified if `learn = "ae"`, a list containing the following named elements to define the architecture and training: `layer_1_dim`, `link_fun`, `epochs`, `loss` and `batch_size`.
#' @param tolerance_level A (typically small) value that we want a quantile (defined by `cutoff_criterion`) of our individual cross-validated losses to be less than (called $\epsilon$ in paper). For example, `tolerance_level = 0.05` and `cutoff_criterion = 0.95` means that we would want 95\% of individual cross-validated losses to be below 0.05. Defaults to 0.05.
#' @param cutoff_criterion A (typically large) quantile (called $\alpha$ in paper) such that we want this quantile of individual cross-validated losses to be less than the value defined by `tolerance_level`. For example, `tolerance_level = 0.05` and `cutoff_criterion = 0.95` means that we would want 95\% of individual cross-validated losses to be below 0.05. Defaults to 0.05.
#' @param learn_function a function only to be supplied if `learn = user`, a user-defined function that takes arguments `Y` (data matrix) and `k` (latent dimension) and returns a list containing two elements `Encode` and `Decode` which define the custom encoding and decoding transformations, respectively.
#' @param verbose logical, whether to print output to console during training and cross-validation. Defaults to `TRUE`.
#'
#' @return
#' @returns
#' @export
#' @importFrom magrittr "%>%"
#'
#' @examples
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
    qc <- NA
  } else {
    index <- min(which(cutoff_criterion_quantiles <= tolerance_level))
    qc <- breaks[index]
  }

  # Loss-lessness (Scree) Plot:  ---------------------------------------------
  summary_correlation_plot(out_basisel = out,
                           cvqlines = cvqlines,
                           cutoff_criterion = cutoff_criterion,
                           r = r,
                           q = q,
                           breaks = breaks,
                           method_name = method_name,
                           qc = qc,
                           tolerance_level = tolerance_level)


  # Heatmap: ----------------------------------------------------------------
 heat_map <- plotly::plot_ly(
    x = breaks,
    y = seq(0, 1, length.out = n),
    z = log10(out$Qrho_v),
    type = "heatmap",
    colors = "RdYlGn",
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
  if(!is.na(qc)) {
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
        learn_final[["Encode"]](Y = Y, k = qc)
      }
    } else {
      learn_final <- learn_function(Y = mat, k = qc)
      out[["Encode"]] <- learn_final[["Encode"]]
    }

    out[["Decode"]] <- learn_final[["Decode"]]

  } else if(is.na(qc)){
    out[["Encode"]] <- NA
    out[["Decode"]] <- NA
  }

  out$qc <- qc
  out$heatmap <- heat_map
  # final:
  out
}
