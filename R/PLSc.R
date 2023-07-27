#' @title PLS
#' @description Performs PLS two class
#' @usage PLSc(X, Y, A, scaling, post.transformation, eps, Y.prob)
#' @param X data matrix where columns represent the \eqn{p} classes and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and rows the \eqn{n} observations.
#' @param A number of components
#' @param scaling type of scaling.
#' @param post.transformation TRUE if you want to apply post transformation.
#' @param eps parameter needed to transform dummy variable into a
#' vector of probabilities
#' @param Y.prob TRUE if Y describes the probability to being in the class.
#' @author Angela Andreella
#' @return Returns a list with the following objects:
#' - \code{W}: matrix of weigths
#' - \code{X_loading}: matrix of X loading
#' - \code{Y_loading}: matirx of Y loading
#' - \code{X}: matrix of X data
#' - \code{Y}: matirx of Y data
#' - \code{T_score}: matrix of scores
#' - \code{Y_fitted}: fitted Y matrix
#' - \code{B}: Matrix regression coefficients
#' - \code{M}: number of orthogonal components if post transformation is applied.
#' @importFrom compositions ilrInv
#' @importFrom compositions ilr
#' @importFrom stats sd
#' @importFrom stats model.matrix
#'
#' @export



PLSc <- function(X, Y, A, scaling, post.transformation,
                 eps, Y.prob){

  nY <- ifelse(is.null(dim(Y)), length(Y), dim(Y)[1])
  if(dim(X)[1] != nY){
    stop("X and Y must have the same number of observations!")
    }

  if(!(scaling %in% c("auto-scaling", "pareto-scaling", "mean-centering"))){
    stop("available scaling are auto-scaling, pareto-scaling and mean-centering")
  }
  X<-as.matrix(X)
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

  #If Y is not probability but a vector of classes
  if(!Y.prob){

    if(is.null(dim(Y)) | ncol(as.matrix(Y))==1){
      Y <- as.matrix(Y)

      if(!is.factor(Y)){
        Y <- as.factor(Y)

      }
      levels(Y) <- c(0,1)
      Y <- model.matrix(~0+Y)
    }

    #Transform to probability matrix
    Y[which(Y==0)]<-eps
    Y[which(Y==1)]<-1-(ncol(Y)-1)*eps

    #Centered log ratio transform transformation
    P <- matrix(clr(Y), ncol = ncol(Y))
  }else{
    P <- Y
  }


  #scaling Y
  Mm <- apply(P, 2, mean)
  s <- apply(P, 2, sd)
  P <- (P - Mm)/s



  n <- nrow(X)

  P <- as.matrix(P)
  X <- as.matrix(X)
  out <- computeWT(X = X, Y = P, A = A)

  W <- out$W
  T_score <- out$T_score

  R <- out$R
  if(post.transformation){

    out <- ptPLSc(X = X, Y = matrix(clr(Y), ncol = ncol(Y)), W = W)
    Wtilde <- out$Wtilde
    M <- out$M

    #apply G to weight matrix
    N <- dim(X)[1]
    E <- r <- Q <- list()

    E[[1]] <- X

    T_score <- IDA(X = X, Y = Y, W = Wtilde)
  }else{
    M <- NULL
  }

  #Compute loadings matrix
  X_loading = t(X) %*% T_score %*% solve(t(T_score) %*% T_score)
  Y_loading = t(Y) %*% T_score %*% solve(t(T_score) %*% T_score)

  #matrix coefficients
  Wstar = W %*% solve(t(X_loading) %*%W)
  B = Wstar %*% t(Y_loading)


  Y_fitted <- fitY(X = X, B = B, Mm = Mm, s = s)



  return(list(X_loading = X_loading,
              Y_loading = Y_loading,
              X = X,
              Y = Y,
              B = B,
              M = M,
              T_score = T_score,
              Y_fitted = Y_fitted))
}
