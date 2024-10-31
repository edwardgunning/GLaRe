library(GLarE)
mnist <- keras::dataset_mnist()
## normalize so the range is (0,1)
x_train <- mnist$train$x/255

tst_learn <- learn_dwt.2d(Y = x_train[1:100,,])
D_thresh <-
Yhat <- tst_learn$Decode(Ystar = tst_learn$Encode(Y = x_train[1:100,, ], k = 32 * ))


par(mfrow = c(1, 2))
image(x_train[10,,])
image(tst_learn$Decode(Ystar = tst_learn$Encode(Y = x_train[1:100,, ], k = 2))[10,, ])



dwt.2d_GLaRE <- GLaRe(mat = x_train[1:100,,],
                      lim = 100,
                      incr = 10,
                      learn = "dwt.2d",
                      kf = 5,
                      sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
                      cvqlines = 0.5,
                      cutoffcriterion = 0.05,
                      cutoffvalue = 0.95,
                      verbose = TRUE)


options(error = recover)

par(mfrow = c(1, 2), cex = 0.5)
pca_glare <- GLaRe(
  mat = x_train,
  learn = "pca",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.05,
  cutoffvalue = 0.95,
  incr = 10,
  lim = 150,
  verbose = TRUE)

ae_glare <- GLaRe(
  mat = x_train,
  learn = "autoencoder",
  kf = 5,
  sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
  cvqlines = 0.5,
  cutoffcriterion = 0.05,
  cutoffvalue = 0.95,
  incr = 10,
  lim = 150,
  ae_args = list(link_fun = "sigmoid", epoch = 50, loss = "binary_crossentropy"))

save.image(file = "demos/mnist.rda")

