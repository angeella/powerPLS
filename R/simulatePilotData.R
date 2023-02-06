#' @title simulate cluster data from normal distribution
#' @description simulate simple cluster data from normal distribution
#' @usage simulatePilotData(seed=123, nvar, clus.size,m, s)
#' @param seed seed value
#' @param nvar number of variables
#' @param clus.size size of classes (only two classes are considered)
#' @param m mean multivariate normal distribution
#' @param s sd multivariate normal distribution
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @author Angela Andreella
#' @return Returns list of X and Y simulated data
#' @export


simulatePilotData <- function(seed=123, nvar, clus.size,m, s){

  #Simulate data as PT + E where P is composed by orthogonal components,
  #T from matrix normal distribution and E noise

  Tm <- rbind(rmvnorm(n = clus.size[1], mean = m),
              rmvnorm(n = clus.size[2], mean = rep(0,nvar)))

  P <- svd(matrix(rnorm((clus.size[1]+clus.size[2])*(clus.size[1]+clus.size[2])),ncol = clus.size[1]+clus.size[2]))$u

  X<- P %*% Tm + matrix(rnorm(nvar*(clus.size[1]+clus.size[2])),ncol = nvar)

  Y<- c(rep(0, clus.size[1]), rep(1, clus.size[2]))

  return(simData = list(X = X, Y = Y))
}
