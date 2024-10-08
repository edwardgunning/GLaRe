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
  cutoffcriterion = 0.05,
  cutoffvalue = 0.95,
  incr = 8,
  lim = 93)

DTI_ae <- GLaRe(
  mat = DTI,
  learn = "autoencoder",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.05,
  cutoffvalue = 0.95,
  incr = 8,
  lim = 93,
  ae_args = list(link_fun = "linear", epoch = 50))

DTI_dwt <- GLaRe(
  mat = DTI,
  learn = "dwt",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.05,
  cutoffvalue = 0.95,
  incr = 8,
  lim = 93)



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
#   lim = 10
# )

PH <- readr::read_table(file = "https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/npfda-phoneme.dat", col_names = FALSE)
PH <- as.matrix(PH)[, 1:150]
matplot(t(PH), type = "l")
par(mfrow = c(1, 2))
test_glare <- GLaRe(
  mat = PH,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.5,
  cutoffvalue = 0.9,
  incr = 5,
  lim = 50)

test_glare <- GLaRe(
  mat = DTI,
  learn = "autoencoder",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.5,
  cutoffvalue = 0.9,
  incr = 5,
  lim = 50,
  ae_args = list(link_fun = "linear", epoch = 50))




