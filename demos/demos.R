library(GLarE)

# Set Random Seed: --------------------------------------------------------
tensorflow::set_random_seed(seed = 1996)
data("glaucoma_data")

glaucoma <- as.matrix(glaucoma_data)
plot_eye(y = glaucoma[1, ])

glaucoma_pca <- GLaRe(
  mat = glaucoma,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 20,
  latent_dim_to = 500)


plot_eye_reconstruction(GLaRe_output = glaucoma_pca, y = glaucoma[50,])


glaucoma_reshaped <- array(glaucoma, dim = c(nrow(glaucoma), 120, 120))
glaucoma_dwt.2d <- GLaRe(
  mat = glaucoma_reshaped,
  learn = "dwt.2d",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 20,
  latent_dim_to = 500)

glaucoma_ae <- GLaRe(
  mat = glaucoma,
  learn = "ae",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 1,
  latent_dim_to = 1,
  ae_args = list(link_fun = "linear", epoch = 50))


# DTI: --------------------------------------------------------------------

DTI <- refund::DTI$cca
DTI <- na.omit(DTI)


par(mfrow = c(1, 3))
par(mfrow = c(1, 2), cex = 0.55)
DTI_pca <- GLaRe(
  mat = DTI,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 2,
  latent_dim_to = 93)

DTI_ae <- GLaRe(
  mat = DTI,
  learn = "ae",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 8,
  latent_dim_to = 93,
  ae_args = list(link_fun = "linear", epoch = 50))

DTI_dwt <- GLaRe(
  mat = DTI,
  learn = "dwt",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 2,
  latent_dim_to = 93)


plot_1D_reconstruction(GLaRe_output = DTI_pca, Y = DTI[1:5,])

par(mfrow = c(1, 2))
plot_train_validation_ratio(GLaRe_output = DTI_dwt)
plot_train_validation_ratio(GLaRe_output = DTI_pca)

par(mfrow = c(1, 2))
distribution_plot(GLaRe_output = DTI_pca)
distribution_plot(GLaRe_output = DTI_dwt)



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
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  verbose = TRUE)

PH_dwt <- GLaRe(
  mat = PH,
  learn = "dwt",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.9,
  cutoff_criterion = 0.95,
  tolerance_level = 0.05,
  latent_dim_by = 10,
  latent_dim_to = ncol(PH),
  verbose = TRUE)



par(mfrow = c(1, 2))
GLarE:::plot_train_validation_ratio(GLaRe_output = PH_dwt)
GLarE:::plot_train_validation_ratio(GLaRe_output = PH_pca)


