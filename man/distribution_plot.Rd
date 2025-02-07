% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper-functions.R
\name{distribution_plot}
\alias{distribution_plot}
\title{Plot Distribution of Validation Losses}
\usage{
distribution_plot(GLaRe_output, jitter = 0.2, cex = 1)
}
\arguments{
\item{GLaRe_output}{Output object from the `GLaRe()` function.}

\item{jitter}{Amount of jitter to add to the points. Defaults to 0.2.}

\item{cex}{Point size for the plot. Defaults to 1.}
}
\value{
A strip chart of validation losses.
}
\description{
Creates a jittered strip chart to visualize the distribution of individual validation losses
for different numbers of latent features.
}
