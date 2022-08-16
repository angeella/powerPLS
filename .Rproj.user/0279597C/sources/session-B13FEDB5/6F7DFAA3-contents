test_that("ptPLSc works", {

  data("umorAcqueo")

  X <- umorAcqueo$X
  Y <- umorAcqueo$Y


  out <- computeWT(X = X, Y = Y, A = 4)
  out <- ptPLSc(X = X, Y = Y, W = out$W)
  M <- out$M
  check1 <- ifelse(is.character(all.equal(round(out$G %*% t(out$G)),diag(4))), FALSE, TRUE)

  out <- PLSc(X = X, Y = Y, A = 4, scaling = "mean-centering", post.transformation = TRUE, eps = 0.01)

  Y[which(Y==0)]<-0.01

  Y[which(Y==1)]<-1-(ncol(Y)-1)*0.01
  P <- matrix(compositions::clr(Y), ncol = ncol(Y))
  check2 <- ifelse(is.character(all.equal(round( t(P) %*% out$T.score[,1:M]), diag(M) - diag(M))), FALSE, TRUE)


  expect_equal(TRUE, check1 & check2)
})
