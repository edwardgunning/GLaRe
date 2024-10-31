transform_correlation_output <- function(out_basissel, cvqlines) {
  if (!(cvqlines >= 0 & cvqlines <= 1)) stop("cvqlines must be in [0, 1]")
  if (!(is.list(out_basissel))) stop("out_basissel must be the named list output")
  if (!all(c("corM_t", "rho_v") %in% names(out_basissel))) stop("out_basissel object must contain elements corM_t and rho_v")
  cor_df <- data.frame(
    meansqcor_t = out_basissel[["corM_t"]],
    minsqcor_cv = apply(out_basissel[["rho_v"]], 2, min),
    meansqcor_cv = apply(out_basissel[["rho_v"]], 2, function(x) sum(x) / (length(x))),
    maxsqcor_cv = apply(out_basissel[["rho_v"]], 2, max),
    medsqcor_cv = apply(out_basissel[["rho_v"]], 2, function(x) quantile(x, 0.5, na.rm = TRUE)),
    qchoice_cv = apply(out_basissel[["rho_v"]], 2, function(x) quantile(x, cvqlines, na.rm = TRUE))
  )
  cor_df
}

summary_correlation_plot <- function(out_basisel, cvqlines, r, q, breaks, method_name, qc, cutoffvalue) {
  correlation_df <- transform_correlation_output(out_basisel, cvqlines)
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

  if (!is.na(qc)) {
    abline(v = qc, lty = 2, col = "grey")
    abline(h = cutoffvalue, lty = 2, col = "grey")
    axis(side = 1, at = c(qc), labels = paste0("qc = ", qc), col = "darkgrey", font = 4, lwd = 3, padj = 1.2)
  }


  legend("topright",
    legend = c(
      bquote(paste("CV Min ", R^2, sep = "")),
      bquote(paste("CV Mean ", R^2, sep = "")),
      bquote(paste("CV Percentile = ", .(cvqlines), " ", R^2, sep = "")),
      bquote(paste("CV Max ", R^2, sep = "")),
      bquote(paste("Train Mean ", R^2, sep = ""))
    ),
    col = c("blue", "goldenrod", "purple", "red3", "green"),
    lty = c(1, 1, 1, 1, 1),
    lwd = c(2, 2, 2, 2, 2, 2)
  )
}




#' Assess losslessness of a latent feature representation tool.
#'
#' @param mat An n-by-p data matrix (n observations, p variables).
#' @param learn The latent feature representation chosen, defaults to principal component analysis (PCA).
#' @param kf The number of folds for the k-fold cross-validation.
#' @param center Whether or not the data should be mean-centered prior to performing the latent feature transformation.
#' @param scale Whether the data should be scaled prior to analysis, so that each variable has standard deviation 1.
#' @param sqcorrel Quantiles of the training and test squared correlation to return.
#' @param cvqlines The user-specified quantile of the cross-validated squared correlation to display on the plot.
#' @param cutoffcriterion
#' @param cutoffvalue
#'
#' @return
#' @export
#' @importFrom magrittr "%>%"
#'
#' @examples
GLaRe <- function(
    mat = iris[, 1:4],
    latent_dim_from = 1,
    latent_dim_to = min(ncol(mat) - 1, nrow(mat) - 1),
    latent_dim_by = 1,
    learn = "pca",
    method_name = toupper(learn),
    kf = 5,
    sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
    cvqlines = .5,
    ae_args = list(),
    cutoffcriterion = 0.5,
    cutoffvalue = 0.9,
    verbose = TRUE) {
  print(paste("*** Learning Method:", learn, "***"))
  if (learn == "pca") {
    out <- flf_basissel_pca(mat = mat, kf = kf, latent_dim_from = 1, latent_dim_to = latent_dim_to, latent_dim_by = latent_dim_by, verbose = verbose)
  } else if (learn == "ae") {
    out <- flf_basissel_ae(mat = mat, kf = kf, latent_dim_from = 1, latent_dim_to = latent_dim_to, latent_dim_by = latent_dim_by, ae_args = ae_args, verbose = verbose)
  } else if (learn == "dwt") {
    out <- flf_basissel_dwt(mat = mat, kf = kf, latent_dim_from = 1, latent_dim_to = latent_dim_to, latent_dim_by = latent_dim_by, verbose = verbose)
  } else if (learn == "dwt.2d") {
    out <- flf_basissel_dwt.2d(mat = mat, kf = kf, latent_dim_from = 1, latent_dim_to = latent_dim_to, latent_dim_by = latent_dim_by, verbose = verbose)
  }

  # else {
  #  #   out <- flf_basissel_user(mat, learn_user, kf, center, scale, ...)
  #  # }
  n <- nrow(mat)
  p <- ncol(mat)
  q <- out[["q"]]
  r <- out[["r"]]
  breaks <- out[["breaks"]]


  # Qualifying Criterion and Add to plot: -----------------------------------
  cutoff_criterion_quantiles <- apply(out[["rho_v"]][, seq_len(r)], 2, function(x) quantile(x, cutoffcriterion, na.rm = TRUE))
  if (!any(cutoff_criterion_quantiles <= cutoffvalue)) {
    warning("No qualifying criterion found, try adjusting parameters.")
    qc <- NA
  } else {
    index <- min(which(cutoff_criterion_quantiles <= cutoffvalue))
    qc <- breaks[index]
  }

  # Loss-lessness (Scree) Plot:  ---------------------------------------------
  summary_correlation_plot(out_basisel = out, cvqlines = cvqlines, r = r, q = q, breaks = breaks, method_name = method_name, qc = qc, cutoffvalue = cutoffvalue)


  # Heatmap: ----------------------------------------------------------------
  p2 <- plotly::plot_ly(
    x = breaks,
    y = seq(0, 1, length.out = n),
    z = out$Qrho_v,
    type = "heatmap", colors = "RdYlGn",
  ) %>%
    plotly::layout(
      title = "Squared Correlation Heatmap",
      xaxis = list(title = "No. of Latent Features"),
      yaxis = list(title = "Quantile of Observations"),
      legend = list(orientation = "h")
    )

  out$qc <- qc
  out$heatmap <- p2
  # final
  out
}
