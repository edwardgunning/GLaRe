library(GLarE)
mnist <- keras::dataset_mnist()
## normalize so the range is (0,1)
x_train <- mnist$train$x/255
x_train <- keras::array_reshape(x_train, c(nrow(x_train), 28*28), order = "F")

par(mfrow = c(1, 2), cex = 0.5)
pca_glare <- GLaRe(
  mat = x_train,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.5,
  cutoffvalue = 0.9,
  incr = 10,
  lim = 51,
  verbose = TRUE)

ae_glare <- GLaRe(
  mat = x_train,
  learn = "autoencoder",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.5,
  cutoffvalue = 0.9,
  incr = 10,
  lim = 51,
  ae_args = list(link_fun = "sigmoid", epoch = 50, loss = "binary_crossentropy"))

save.image(file = "demos/mnist.rda")

