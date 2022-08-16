test_that("eigenTest works", {

  data("umorAcqueo")

  X <- umorAcqueo$X
  Y <- umorAcqueo$Y

  out <- eigenTest(X = X, Y = Y, nperm = 1000, A = 5, Y.prob = FALSE, X.scaling = FALSE)

  out1 <- mdatools::plsda(x = X, c = as.factor(Y[,1]), ncomp = 4, cv =5)

  check1 <- sum(sort(out$test, decreasing = TRUE) == out$test) == length(out$test)
  check2 <- sum(out$pv_adj <=0.05) == out1$ncomp.selected

  expect_equal(check1 & check2, TRUE)
})
