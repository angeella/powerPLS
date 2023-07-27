#' @title mcc test
#' @description Performs randomization test based on MCC
#' @usage mccTest(X, Y, nperm, A, randomization, ...)
#' @param X data matrix where columns represent the \eqn{p} classes and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and
#' rows the \eqn{n} observations.
#' @param nperm number of permutations
#' @param A number of latent components
#' @param randomization compute p-values?
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns a list with the corresponding statistical tests,
#' raw and adjusted p-values
#' @importFrom compositions clrInv
#' @importFrom compositions clr
#' @importFrom stats sd
#' @export




mccTest <- function(X, Y, nperm = 100, A, randomization = FALSE, ...){

  out <- PLSc(X = X, Y = Y, A = A, ...)
  if(!is.null(dim(out$Y_fitted))){
    Y_fitted <- as.factor(out$Y_fitted[,2])
  }else{
    Y_fitted <- as.factor(out$Y_fitted)
  }
  Yf <- as.factor(Y)
  levels(Yf) <-  c(0,1)
  levels(Y_fitted) <-  c(0,1)
  confMatrix <- table(Yf, Y_fitted)
  mcc_obs <- mcc(confMatrix = confMatrix)

  if(randomization){
    pv <- 0
    for(j in seq(nperm-1)){
      idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
      Xkp <- X[idx,]

      out <- PLSc(X = Xkp, Y = Y, A = A, ...)

      if(!is.null(dim(out$Y_fitted))){
        Y_fitted <- as.factor(out$Y_fitted[,2])
      }else{
        Y_fitted <- as.factor(out$Y_fitted)
      }
      Yf <- as.factor(Y)
      levels(Yf) <- c(0,1)
      levels(Y_fitted) <- c(0,1)
      confMatrix <- table(Yf, Y_fitted)
      mcc_p <- mcc(confMatrix = confMatrix)

      if(mcc_p >= mcc_obs){pv <- pv+1}

    }

    pv <- (pv+1)/nperm
  }else{
    pv <- NA
  }


  return(list(pv = pv,
              pv_adj = min(pv*A,1), test = mcc_obs))
}
