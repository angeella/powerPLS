#' @export
#' @importFrom mvtnorm rmvnorm
#'


sim_pilotData <- function(n, p, q, A, structured.noise, d.class = 10, prop = 0.5){


  k <- A+p
  n1 <- round(n*prop)
  n2 <- n - n1
  #Simulo la matrice degli scores
  rand_matA <- mvtnorm::rmvnorm(n1, mean = rep(d.class,k),
                                      sigma = diag(k))

  rand_matB <- mvtnorm::rmvnorm(n2, mean = rep(-d.class,k), sigma = diag(k))

  rand_mat <- rbind(rand_matA, rand_matB) #X matrix

  out <- prcomp(rand_mat)

  Tprime <- rand_mat %*% out$rotation

  Tscore <- Tprime[,1:A]

  Tscore[,A] <- Tscore[,A]

  rand_mat <- matrix(runif(A*(q+p)), nrow = A)
  out1 <- prcomp(rand_mat)

  X <- Tscore %*% t(out1$rotation)

  X[,1:p] <- X[,1:p] + Tprime[,(A+1):k]*0.1

  Y <- X[,1:p]
  X <- X[,(p+1):q]

  return(list(X = X, Y = Y))
}
