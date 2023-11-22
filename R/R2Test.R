#' @title R2 test
#' @description Performs randomization test based on R2
#' @usage R2Test(X, Y, nperm = 100, A, randomization = FALSE, Y.prob = FALSE, eps,...)
#' @param X data matrix where columns represent the \eqn{p} classes and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and
#' rows the \eqn{n} observations.
#' @param nperm number of permutations
#' @param A number of latent components
#' @param randomization compute p-values?
#' @param Y.prob Is Y a probability vector?
#' @param eps parameter needed to transform dummy variable into a
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns a list with the corresponding statistical tests,
#' raw and adjusted p-values
#' @importFrom compositions ilr
#' @importFrom stats cor
#' @export

R2Test <- function(X, Y, nperm = 100, A, randomization = FALSE, Y.prob = FALSE, eps,...){


  out <- PLSc(X = X, Y = Y, A = A, Y.prob = Y.prob, transformation = "ilr", eps = eps, ...)

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
    P <- matrix(ilr(Y), ncol = 1)
  }else{
    P <- Y
  }
  s <- sd(P)
  Yfitted = matrix(ilrInv(s*(X %*% out$B))[,3], ncol = 1)

  r2_obs <- cor(Yfitted, P)

  if(randomization){
    pv <- 0
    for(j in seq(nperm-1)){
      idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
      Xkp <- X[idx,]

      out <- PLSc(X = Xkp, Y = Y, A = A, Y.prob = Y.prob, eps = eps,...)

      Yfitted =matrix(ilrInv(s*(Xkp %*% out$B))[,3], ncol = 1)


      r2_p <- cor(Yfitted, P)


      if(r2_p >= r2_obs){pv <- pv+1}

    }

    pv <- (pv+1)/nperm
  }else{
    pv <- NA
  }


  return(list(pv = pv,
              pv_adj = min(pv*A,1), test = r2_obs))

}
