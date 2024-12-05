library(GLarE)
mnist <- keras::dataset_mnist()
## normalize so the range is (0,1)
x_train <- mnist$train$x/255
# n <- nrow(x_train)
# x_train <- x_train[1:n,,]
x_train_flattened <- matrix(x_train, nrow(x_train), 784)
mnist_pca <- readRDS( "demos/mnist_pca.rds")
mnist_pca$qc
plot_mnist_reconstruction(GLaRe_output = mnist_pca, y = x_train_flattened[2,])

par(mfrow = c(1, 2), cex = 0.5)

mnist_dwt.2d <- GLaRe(mat = x_train,
                      latent_dim_from = 1,
                      latent_dim_to = 400,
                      latent_dim_by = 20,
                      learn = "dwt.2d",
                      kf = 5,
                      sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
                      cvqlines = 0.5,
                      cutoff_criterion = 0.95,
                      tolerance_level = 0.05,
                      verbose = TRUE)

mnist_pca <- GLaRe(mat = x_train_flattened,
                latent_dim_from = 1,
                latent_dim_to = 1000,
                latent_dim_by = 100,
                learn = "pca",
                kf = 5,
                sqcorrel = c("trainmean", "cvmean", "cvmin", "cvmax"),
                cvqlines = 0.5,
                cutoff_criterion = 0.95,
                tolerance_level = 0.05,
                verbose = TRUE)
saveRDS(object = mnist_pca, file = "demos/mnist_pca.rds")



