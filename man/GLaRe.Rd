% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GlaRe.R
\name{GLaRe}
\alias{GLaRe}
\title{Assess losslessness of a latent feature representation tool.}
\usage{
GLaRe(
  mat = iris[, 1:4],
  lim = min(ncol(mat) - 1, nrow(mat) - 1),
  incr = 1,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  ae_args = list(),
  cutoffcriterion = 0.5,
  cutoffvalue = 0.9,
  verbose = TRUE
)
}
\arguments{
\item{mat}{An n-by-p data matrix (n observations, p variables).}

\item{learn}{The latent feature representation chosen, defaults to principal component analysis (PCA).}

\item{kf}{The number of folds for the k-fold cross-validation.}

\item{sqcorrel}{Quantiles of the training and test squared correlation to return.}

\item{cvqlines}{The user-specified quantile of the cross-validated squared correlation to display on the plot.}

\item{cutoffvalue}{}

\item{center}{Whether or not the data should be mean-centered prior to performing the latent feature transformation.}

\item{scale}{Whether the data should be scaled prior to analysis, so that each variable has standard deviation 1.}
}
\description{
Assess losslessness of a latent feature representation tool.
}
