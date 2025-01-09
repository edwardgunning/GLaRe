
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
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  verbose = FALSE
)

## autoencoder
# ph_ae <- GLaRe(
#   mat = PH,
#   learn = "ae",
#   kf = 5,
#   cvqlines = 0.9,
#   cutoff_criterion = 0.95,
#   tolerance_level = 0.05,
#   latent_dim_by = 10,
#   latent_dim_to = ncol(PH),
#   ae_args = list(link_fun = "linear", epoch = 50),
#   verbose = FALSE
# )

## dwt
ph_dwt <- GLaRe(
  mat = PH,
  learn = "dwt",
  kf = 5,
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  verbose = FALSE
)
#> Warning in GLaRe(mat = PH, learn = "dwt", kf = 5, cvqlines = 0.9,
#> cutoff_criterion = 0.95, : No qualifying criterion found, try adjusting
#> parameters.
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Look at individual losses distributions:

``` r
par(mfrow = c(1, 3))
distribution_plot(GLaRe_output = ph_pca)
# distribution_plot(GLaRe_output = ph_ae)
distribution_plot(GLaRe_output = ph_dwt)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Look at the training validation ratios:

``` r
par(mfrow = c(1, 3))
plot_train_validation_ratio(GLaRe_output = ph_pca)
# plot_train_validation_ratio(GLaRe_output = ph_ae)
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
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 8,
  latent_dim_to = nrow(DTI),
  verbose = FALSE
)

# DTI_ae <- GLaRe(
#   mat = DTI,
#   learn = "ae",
#   kf = 5,
#   cvqlines = 0.9,
#   cutoff_criterion = 0.95,
#   tolerance_level = 0.05,
#   latent_dim_by = 8,
#   latent_dim_to = nrow(DTI),
#   ae_args = list(link_fun = "linear", epoch = 50),
#   verbose = FALSE
# )

DTI_dwt <- GLaRe(
  mat = DTI,
  learn = "dwt",
  kf = 5,
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 8,
  latent_dim_to = nrow(DTI),
  verbose = FALSE
)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

Look at individual losses distributions:

``` r
par(mfrow = c(1, 3))
distribution_plot(GLaRe_output = DTI_pca)
# distribution_plot(GLaRe_output = DTI_ae)
distribution_plot(GLaRe_output = DTI_dwt)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

Look at the training validation ratios:

``` r
par(mfrow = c(1, 3))
plot_train_validation_ratio(GLaRe_output = DTI_pca)
# plot_train_validation_ratio(GLaRe_output = DTI_ae)
plot_train_validation_ratio(GLaRe_output = DTI_dwt)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

## Example: The Glaucoma Data

This data is stored in the package as `glaucoma_data`, along with a
function to plot it, known as `plot_eye()`:

``` r
glaucoma <- as.matrix(glaucoma_data)
plot_eye(y = glaucoma[1, ])
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

``` r
par(mfrow = c(1, 3), cex = 0.5) # side-by-side plots
glaucoma_pca <- GLaRe(
  mat = glaucoma,
  learn = "pca",
  kf = 5,
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 20,
  latent_dim_to = 500)
#> [1] "*** Learning Method: pca ***"
#> [1] "Number of non-zero eigenvectors is less than latent_dim_to"
#> [1] "====== Training ======"
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "= Latent Dim. = 261"
#> [1] "====== Performing 5-fold CV ======="
#> [1] "==== Fold  1 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "==== Fold  2 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "==== Fold  3 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "==== Fold  4 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "==== Fold  5 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "====== Finished 5-fold CV, Summarising Results ======="
#> [1] "Final training of Model at qualifying criterion:"
```

``` r

glaucoma_reshaped <- array(glaucoma, dim = c(nrow(glaucoma), 120, 120))
glaucoma_dwt.2d <- GLaRe(
  mat = glaucoma_reshaped,
  learn = "dwt.2d",
  kf = 5,
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 20,
  latent_dim_to = 500)
#> [1] "*** Learning Method: dwt.2d ***"
#> [1] "====== Training ======"
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "= Latent Dim. = 261"
#> [1] "= Latent Dim. = 281"
#> [1] "= Latent Dim. = 301"
#> [1] "= Latent Dim. = 321"
#> [1] "= Latent Dim. = 341"
#> [1] "= Latent Dim. = 361"
#> [1] "= Latent Dim. = 381"
#> [1] "= Latent Dim. = 401"
#> [1] "= Latent Dim. = 421"
#> [1] "= Latent Dim. = 441"
#> [1] "= Latent Dim. = 461"
#> [1] "= Latent Dim. = 481"
#> [1] "====== Performing 5-fold CV ======="
#> [1] "==== Fold  1 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "= Latent Dim. = 261"
#> [1] "= Latent Dim. = 281"
#> [1] "= Latent Dim. = 301"
#> [1] "= Latent Dim. = 321"
#> [1] "= Latent Dim. = 341"
#> [1] "= Latent Dim. = 361"
#> [1] "= Latent Dim. = 381"
#> [1] "= Latent Dim. = 401"
#> [1] "= Latent Dim. = 421"
#> [1] "= Latent Dim. = 441"
#> [1] "= Latent Dim. = 461"
#> [1] "= Latent Dim. = 481"
#> [1] "==== Fold  2 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "= Latent Dim. = 261"
#> [1] "= Latent Dim. = 281"
#> [1] "= Latent Dim. = 301"
#> [1] "= Latent Dim. = 321"
#> [1] "= Latent Dim. = 341"
#> [1] "= Latent Dim. = 361"
#> [1] "= Latent Dim. = 381"
#> [1] "= Latent Dim. = 401"
#> [1] "= Latent Dim. = 421"
#> [1] "= Latent Dim. = 441"
#> [1] "= Latent Dim. = 461"
#> [1] "= Latent Dim. = 481"
#> [1] "==== Fold  3 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "= Latent Dim. = 261"
#> [1] "= Latent Dim. = 281"
#> [1] "= Latent Dim. = 301"
#> [1] "= Latent Dim. = 321"
#> [1] "= Latent Dim. = 341"
#> [1] "= Latent Dim. = 361"
#> [1] "= Latent Dim. = 381"
#> [1] "= Latent Dim. = 401"
#> [1] "= Latent Dim. = 421"
#> [1] "= Latent Dim. = 441"
#> [1] "= Latent Dim. = 461"
#> [1] "= Latent Dim. = 481"
#> [1] "==== Fold  4 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "= Latent Dim. = 261"
#> [1] "= Latent Dim. = 281"
#> [1] "= Latent Dim. = 301"
#> [1] "= Latent Dim. = 321"
#> [1] "= Latent Dim. = 341"
#> [1] "= Latent Dim. = 361"
#> [1] "= Latent Dim. = 381"
#> [1] "= Latent Dim. = 401"
#> [1] "= Latent Dim. = 421"
#> [1] "= Latent Dim. = 441"
#> [1] "= Latent Dim. = 461"
#> [1] "= Latent Dim. = 481"
#> [1] "==== Fold  5 ===="
#> [1] "= Latent Dim. = 1"
#> [1] "= Latent Dim. = 21"
#> [1] "= Latent Dim. = 41"
#> [1] "= Latent Dim. = 61"
#> [1] "= Latent Dim. = 81"
#> [1] "= Latent Dim. = 101"
#> [1] "= Latent Dim. = 121"
#> [1] "= Latent Dim. = 141"
#> [1] "= Latent Dim. = 161"
#> [1] "= Latent Dim. = 181"
#> [1] "= Latent Dim. = 201"
#> [1] "= Latent Dim. = 221"
#> [1] "= Latent Dim. = 241"
#> [1] "= Latent Dim. = 261"
#> [1] "= Latent Dim. = 281"
#> [1] "= Latent Dim. = 301"
#> [1] "= Latent Dim. = 321"
#> [1] "= Latent Dim. = 341"
#> [1] "= Latent Dim. = 361"
#> [1] "= Latent Dim. = 381"
#> [1] "= Latent Dim. = 401"
#> [1] "= Latent Dim. = 421"
#> [1] "= Latent Dim. = 441"
#> [1] "= Latent Dim. = 461"
#> [1] "= Latent Dim. = 481"
#> [1] "====== Finished 5-fold CV, Summarising Results ======="
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> Warning in loss_function(observed = as.matrix(mat)[i, ], predicted =
#> as.matrix(x)[i, : Predicted values constant: Setting squared correlation to 0
#> [1] "Final training of Model at qualifying criterion:"
```

``` r

# glaucoma_ae <- GLaRe(
#   mat = glaucoma,
#   learn = "ae",
#   kf = 5,
#   cvqlines = 0.9,
#   cutoff_criterion = 0.95,
#   tolerance_level = 0.05,
#   latent_dim_by = 20,
#   latent_dim_to = 500,
#   ae_args = list(link_fun = "linear", epoch = 50))
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

Look at individual losses distributions:

``` r
par(mfrow = c(1, 3))
distribution_plot(GLaRe_output = glaucoma_pca)
# distribution_plot(GLaRe_output = glaucoma_ae)
distribution_plot(GLaRe_output = glaucoma_dwt.2d)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

Look at the training validation ratios:

``` r
par(mfrow = c(1, 3))
plot_train_validation_ratio(GLaRe_output = glaucoma_pca)
# plot_train_validation_ratio(GLaRe_output = glaucoma_ae)
plot_train_validation_ratio(GLaRe_output = glaucoma_dwt.2d)
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

## Example: MNIST

``` r
mnist <- keras::dataset_mnist()
## normalize so the range is (0,1)
mnist_train <- mnist$train$x/255
mnist_train_flattened <- matrix(mnist_train, nrow(mnist_train), 784)
```

``` r
par(mfrow = c(1, 3), cex = 0.5) # side-by-side plots
mnist_pca <- GLaRe(
  mat = mnist_train_flattened,
  learn = "pca",
  kf = 5,
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 20,
  latent_dim_to = 500)

mnist_dwt.2d <- GLaRe(
  mat = mnist_train,
  learn = "dwt.2d",
  kf = 5,
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 20,
  latent_dim_to = 500)

# mnist_ae <- GLaRe(
#   mat = mnist_train_flattened,
#   learn = "ae",
#   kf = 5,
#   cvqlines = 0.9,
#   cutoff_criterion = 0.95,
#   tolerance_level = 0.05,
#   latent_dim_by = 20,
#   latent_dim_to = 500,
#   ae_args = list(link_fun = "linear", epoch = 50))
```

Look at individual losses distributions:

``` r
par(mfrow = c(1, 3))
distribution_plot(GLaRe_output = mnist_pca)
# distribution_plot(GLaRe_output = mnist_ae)
distribution_plot(GLaRe_output = mnist_dwt)
```

Look at the training validation ratios:

``` r
par(mfrow = c(1, 3))
plot_train_validation_ratio(GLaRe_output = mnist_pca)
# plot_train_validation_ratio(GLaRe_output = mnist_ae)
plot_train_validation_ratio(GLaRe_output = mnist_dwt)
```
