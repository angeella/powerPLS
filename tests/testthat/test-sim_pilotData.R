test_that("simulate_pilotData works", {

  out <- sim_pilotData(n = 50, p = 2, q = 60, A = 2, structured.noise = TRUE, d.class = 10, prop = 0.5)

  X <- out$X
  Y <- out$Y

  out <- eigenTest(X = X, Y = Y, nperm = 1000, A = 4, Y.prob = TRUE, X.scaling = "mean-centering")

  check1 <- sum(sort(out$test, decreasing = TRUE)== out$test) == 4


  expect_equal(check1, TRUE)
})
