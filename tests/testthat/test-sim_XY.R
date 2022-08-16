test_that("X-simulation works", {

  data("umorAcqueo")

  X <- umorAcqueo$X
  Y <- umorAcqueo$Y

  out <- PLSc(X = X, Y = Y, A = 4, scaling = "mean-centering", post.transformation = TRUE, eps = 0.01)

  out <- sim_XY(out = out, n = nrow(X), seed = 123, post.transformation = TRUE, A = 4)

  outX <- similarityMatrix(X = X, Y = out$X_H1)
  eigen_d <- ifelse(is.character(all.equal(0, sum(round(outX$eigen_d)), tolerance = 3)), FALSE, TRUE)
  RV <- round(outX$RV)>0.5

  expect_equal(TRUE, eigen_d | RV)
})
