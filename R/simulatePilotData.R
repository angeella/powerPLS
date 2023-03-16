#' @title simulate cluster data from normal distribution
#' @description simulate simple cluster data from normal distribution
#' @usage simulatePilotData(seed=123, nvar, clus.size, m, rho)
#' @param seed seed value
#' @param nvar number of variables
#' @param clus.size size of classes (only two classes are considered)
#' @param m mean multivariate normal distribution
#' @param rho numeric value in `[0,1]`. Level of equi-correlation between pairs
#' of variables
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @author Angela Andreella
#' @return Returns list of X and Y simulated data
#' @export


simulatePilotData <- function(seed=123, nvar, clus.size, m, rho){

  #Simulate data as PT + E where P is composed by orthogonal components,
  #T from matrix normal distribution and E noise
  n <- clus.size[1] + clus.size[2]
  eps <- sqrt(1-rho)*matrix(rnorm(nvar*n), ncol=nvar) + sqrt(rho)*matrix(rep(rnorm(n),nvar), ncol=nvar)
  mu <- c(rep(m,clus.size[1]), rep(0, clus.size[2]))
  Tm <- matrix(rep(mu, nvar), ncol=nvar) + eps

  P <- svd(matrix(rnorm((clus.size[1]+clus.size[2])*(clus.size[1]+clus.size[2])),
                  ncol = clus.size[1]+clus.size[2]))$u

  X<- P %*% Tm

  Y<- c(rep(0, clus.size[1]), rep(1, clus.size[2]))

  return(simData = list(X = X, Y = Y))
}


