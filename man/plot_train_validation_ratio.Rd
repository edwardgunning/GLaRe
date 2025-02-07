% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper-functions.R
\name{plot_train_validation_ratio}
\alias{plot_train_validation_ratio}
\title{Plot Training vs Validation Loss Ratio}
\usage{
plot_train_validation_ratio(GLaRe_output)
}
\arguments{
\item{GLaRe_output}{Output object from the `GLaRe()` function.}
}
\value{
A plot showing the training-to-validation loss ratio.
}
\description{
Visualizes the ratio of average training loss to average validation loss as a function
of the number of latent features.
}
