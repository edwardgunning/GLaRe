library(GLarE)

# Set Random Seed: --------------------------------------------------------
tensorflow::set_random_seed(seed = 1996)

library(fda)

learn_Bspline <- function(Y, k) {
  bspline_basis <- create.bspline.basis(rangeval = c(1, ncol(Y)), nbasis = k, norder = 4)

  Encode <- function(Y) {
    fd_obj <- Data2fd(argvals = 1:ncol(Y), y = t(Y), basisobj = bspline_basis, lambda = 10^-12)
    Ystar <- fd_obj[["coefs"]]
    t(Ystar)
  }

  Decode <- function(Ystar) {
    t(eval.fd(evalarg = 1:ncol(Y), fdobj = fd(coef = t(Ystar), basisobj = bspline_basis)))
  }

  list(Encode = Encode, Decode = Decode)
}

DTI <- refund::DTI$cca
DTI <- na.omit(DTI)

DTI_user <- GLaRe(
  mat = DTI,
  learn = "user",
  kf = 2,
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_from = 4,
  latent_dim_by = 1,
  latent_dim_to = 93,
  learn_function = learn_Bspline,
  method_name = "B-splines")

plot_1D_reconstruction(DT)

DTI_user <- GLaRe(
  mat = DTI,
  learn = "pca",
  kf = 20,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_from = 4,
  latent_dim_by = 1,
  latent_dim_to = 93)

DTI_user <- GLaRe(
  mat = DTI,
  learn = "dwt",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_from = 4,
  latent_dim_by = 1,
  latent_dim_to = 93)



