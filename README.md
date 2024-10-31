
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<center>
<img src="inst/figures/hex.png" alt="hex" width="200"/>
</center>
<!-- badges: end -->

The goal of `GLaRE` is to facilitate the evaluation of losslessness of
different latent feature representations based on generalisation error.

## Installation

You can install the development version of GLarE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("edwardgunning/GLaRe")
```

Load the `GLaRE` package.

``` r
library(GLarE)
```

Set random seed for both shuffling folds and stochastic optimisation of
the autoencoder model.

``` r
tensorflow::set_random_seed(1996)
```

## Example: Phenoeme Data

We use the Phenoeme dataset for this example (see
[here](https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/)).

``` r
# Phoenome dataset:
PH_path <- "https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/npfda-phoneme.dat"
PH <- readr::read_table(file = PH_path, col_names = FALSE)
PH <- as.matrix(PH)[, 1:150]
par(mfrow = c(1, 1))
matplot(t(PH)[, sample(1:nrow(PH), size = 20)], type = "l", xlab = "Freq.", ylab = "Log-periodogram")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="50%" />

We look at the representation given by PCA, autoencoder and DWT.

``` r
# run GLaRe
par(mfrow = c(1, 3), cex = 0.5) # side-by-side plots

## PCA:
ph_pca <- GLaRe(
  mat = PH,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  verbose = FALSE
)
#> [1] "*** Learning Method: pca ***"
```

``` r

## autoencoder
ph_ae <- GLaRe(
  mat = PH,
  learn = "ae",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  ae_args = list(link_fun = "linear", epoch = 50),
  verbose = FALSE
)
#> [1] "*** Learning Method: ae ***"
#> Warning in GLaRe(mat = PH, learn = "ae", kf = 5, sqcorrel = c("trainmean", : No
#> qualifying criterion found, try adjusting parameters.
```

``` r

## dwt
ph_dwt <- GLaRe(
  mat = PH,
  learn = "dwt",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  verbose = FALSE
)
#> [1] "*** Learning Method: dwt ***"
#> Warning in GLaRe(mat = PH, learn = "dwt", kf = 5, sqcorrel = c("trainmean", :
#> No qualifying criterion found, try adjusting parameters.
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Look at individual losses distributions:

``` r
par(mfrow = c(1, 3))
distribution_plot(GLaRe_output = ph_pca)
distribution_plot(GLaRe_output = ph_ae)
distribution_plot(GLaRe_output = ph_dwt)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Look at the training validation ratios:

``` r
par(mfrow = c(1, 3))
plot_train_validation_ratio(GLaRe_output = ph_pca)
plot_train_validation_ratio(GLaRe_output = ph_ae)
plot_train_validation_ratio(GLaRe_output = ph_dwt)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

## Example: DTI Data

From the `refund` R package.

``` r
DTI <- refund::DTI$cca
DTI <- na.omit(DTI)
matplot(t(DTI), type = "l")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="50%" />

``` r
par(mfrow = c(1, 3), cex = 0.5) # side-by-side plots
DTI_pca <- GLaRe(
  mat = DTI,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 8,
  latent_dim_to = nrow(DTI),
  verbose = FALSE
)
#> [1] "*** Learning Method: pca ***"
```

``` r

DTI_ae <- GLaRe(
  mat = DTI,
  learn = "ae",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 8,
  latent_dim_to = nrow(DTI),
  ae_args = list(link_fun = "linear", epoch = 50),
  verbose = FALSE
)
#> [1] "*** Learning Method: ae ***"
#> Warning in GLaRe(mat = DTI, learn = "ae", kf = 5, sqcorrel = c("trainmean", :
#> No qualifying criterion found, try adjusting parameters.
```

``` r

DTI_dwt <- GLaRe(
  mat = DTI,
  learn = "dwt",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 8,
  latent_dim_to = nrow(DTI),
  verbose = FALSE
)
#> [1] "*** Learning Method: dwt ***"
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

Look at individual losses distributions:

``` r
par(mfrow = c(1, 3))
distribution_plot(GLaRe_output = DTI_pca)
distribution_plot(GLaRe_output = DTI_ae)
distribution_plot(GLaRe_output = DTI_dwt)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

Look at the training validation ratios:

``` r
par(mfrow = c(1, 3))
plot_train_validation_ratio(GLaRe_output = DTI_pca)
plot_train_validation_ratio(GLaRe_output = DTI_ae)
plot_train_validation_ratio(GLaRe_output = DTI_dwt)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />
