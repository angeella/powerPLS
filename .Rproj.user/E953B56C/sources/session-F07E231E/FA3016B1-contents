#' @title PLS
#' @description Performs PLS two class
#' @usage PLSc(X, Y, A, scaling, post.transformation = TRUE)
#' @param X data matrix where columns represent the \eqn{p} classes and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and rows the \eqn{n} observations.
#' @param A number of components
#' @author Angela Andreella
#' @return Returns a list with the following objects:
#' - \code{W}: matrix of weigths
#' - \code{X.loading}: matrix of X loading
#' - \code{Y.loading}: matirx of Y loading
#' - \code{T_score}: matrix of scores
#' - \code{Y.fitted}: fitted Y matrix
#' - \code{B}: Matrix regression coefficients
#' - \code{M}: number of orthogonal components if post transformation is applied.
#' @importFrom compositions ilrInv
#' @importFrom compositions ilr
#' @export



PLSc <- function(X, Y, A, scaling, post.transformation = TRUE, eps){

  if(dim(X)[1] != dim(Y)[1]){
    stop("X and Y must have the same number of observations!")
    }

  #TODO check arguments scaling

  if(scaling == "auto-scaling"){
    Mm <- apply(X, 2, mean)
    s <- apply(X, 2, sd)
    X <- (X - Mm)/s
  }
  if(scaling == "pareto-scaling"){
    Mm <- 0
    s <- apply(X, 2, sd)
    X <- (X - Mm)/s
  }
  if(scaling == "mean-centering"){
    Mm <- apply(X, 2, mean)
    s <- 1
    X <- (X - Mm)/s
  }

  Y[which(Y==0)]<-eps

  Y[which(Y==1)]<-1-(ncol(Y)-1)*eps


  #TODO::check if Y must be transformed
  P <- matrix(clr(Y), ncol = ncol(Y))


  Mm <- apply(P, 2, mean)
  s <- apply(P, 2, sd)
  P <- (P - Mm)/s
  #scaling Y


  n <- nrow(X)

  P <- as.matrix(P)
  out <- computeWT(X = X, Y = P, A = A)

  W <- out$W
  T.score <- out$T.score

  R <- out$R
  if(post.transformation){

    out <- ptPLSc(X = X, Y = matrix(clr(Y), ncol = ncol(Y)), W = W)
    Wtilde <- out$Wtilde
    M <- out$M

    #apply G to weight matrix
    N <- dim(X)[1]
    E <- r <- Q <- list()

    E[[1]] <- X

    T.score <- IDA(X = X, Y = Y, W = Wtilde)
  }else{
    M <- NULL
  }

  #Compute loadings matrix
  X.loading = t(X) %*% T.score %*% solve(t(T.score) %*% T.score)
  Y.loading = t(Y) %*% T.score %*% solve(t(T.score) %*% T.score)

  #matrix coefficients
  Wstar = W %*% solve(t(X.loading) %*%W)
  B = Wstar %*% t(Y.loading)

  Y.fitted <- fitY(X = X, B = B, Mm = Mm, s = s)

  return(list(X.loading = X.loading,
              Y.loading = Y.loading,
              B = B,
              M = M,
              T.score = T.score,
              Y.fitted = Y.fitted))
}
