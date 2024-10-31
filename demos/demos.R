library(GLarE)
DTI <- refund::DTI$cca
DTI <- na.omit(DTI)

# Set Random Seed: --------------------------------------------------------
tensorflow::set_random_seed(1996)

par(mfrow = c(1, 3))

DTI_pca <- GLaRe(
  mat = DTI,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 8,
  latent_dim_to = 93)


DTI_ae <- GLaRe(
  mat = DTI,
  learn = "ae",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 8,
  latent_dim_to = 93,
  ae_args = list(link_fun = "linear", epoch = 50))

DTI_dwt <- GLaRe(
  mat = DTI,
  learn = "dwt",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 8,
  latent_dim_to = 93)

par(mfrow = c(1, 2))
plot_train_validation_ratio(GLaRe_output = DTI_dwt)
plot_train_validation_ratio(GLaRe_output = DTI_pca)

par(mfrow = c(1, 2))
distribution_plot(GLaRe_output = DTI_pca)
distribution_plot(GLaRe_output = DTI_dwt)
# Spec <- readr::read_table(file = "https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/npfda-spectrometric.dat", col_names = FALSE)
# Spec <- as.matrix(Spec)[, 1:101]
#
# test_glare_2 <- GLaRe(
#   mat = Spec,
#   learn = "pca",
#   kf = 5,
#   sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
#   cvqlines = 0.5,
#   cutoffcriterion = 3,
#   cutoffvalue = 0.9,
#   latent_dim_to = 10
# )

PH <- readr::read_table(file = "https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/npfda-phoneme.dat", col_names = FALSE)
dim(PH)
PH <- as.matrix(PH)[, 1:150]
matplot(t(PH), type = "l")
par(mfrow = c(1, 2))
PH_pca <- GLaRe(
  mat = PH,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  verbose = TRUE)

PH_dwt <- GLaRe(
  mat = PH,
  learn = "dwt",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.95,
  cutoffvalue = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  verbose = TRUE)




par(mfrow = c(1, 2))
GLarE:::plot_train_validation_ratio(GLaRe_output = PH_dwt)
GLarE:::plot_train_validation_ratio(GLaRe_output = PH_pca)

