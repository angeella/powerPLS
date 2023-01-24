#' @title post transformed PLS
#' @description Performs post transformed PLS
#' @usage ptPLSc(X, Y, W)
#' @param X data matrix where columns represent the \eqn{p} classes and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and rows the \eqn{n} observations.
#' @param W weight matrix where columns represent the \eqn{A} components and rows the \eqn{k} X variables.
#' @author Angela Andreella
#' @return Returns a matrix of scores vectors \code{Tscore}.
#' @export
#' @importFrom compositions clrInv
#' @return Returns a list with the following objects:
#' - \code{W}: matrix of weigths
#' - \code{6}: post transformation matrix
#' - \code{M}: number of orthogonal components.

ptPLSc <- function(X, Y, W){


  #TODO::put check if Y is transformed.

  A <- ncol(W)
  out <- svd(t(Y) %*% X %*% W)
  V <- out$v
  M <- sum(Re(out$d) > 10^-8)
  V <- as.matrix(V[,1:M])
  d2<-eigen((diag(nrow(V))-V%*%t(V))%*%(t(W)%*%t(X)%*%X%*%W))
  #Compute orthogonal part
  Go<-Re(d2$vectors[,1:M])
  lp <- c()
  G<-Go
  for (i in 1:(A-M)) {
    Xp<-X-X%*%W%*%G%*%solve(t(G)%*%t(W)%*%t(X)%*%X%*%W%*%G)%*%t(G)%*%t(W)%*%t(X)%*%X
    d3<-eigen(t(diag(A)-Go%*%t(Go))%*%(t(W)%*%t(Xp)%*%Y%*%t(Y)%*%Xp%*%W))
    Gp<-Re(d3$vectors[,1])
    lp[i]<-d3$values[1]
    G<-cbind(G,Gp)
  }

  if((A-M) <= min(qr(Y)$rank, A)){
    warning("The minimum number of predictive latent variables is greater than min(rank(Y), A)")
  }


  #apply G to weight matrix
  Wtilde <- W %*% G #nX times ncomp

  return(list(Wtilde = Wtilde,
              G = G,
              M = M))

}
