#' @title Optimal number of components by cross-validation
#' @description Performs cross validation to find the optimal number of
#' components by MCC
#' @usage crossValidation(X, Y, K = 10, method = "cv", times = 1,seed = 123,
#' A = 5, randomization = FALSE, ...)
#' @param X data matrix where columns represent the \eqn{p} variables and rows
#' the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} classes and rows
#' the \eqn{n} observations.
#' @param K number of folds
#' @param method \code{"cv"} or \code{"repeated-cv"}
#' @param times if \code{"method = repeated-cv"} here you specify the times
#' @param seed seed value
#' @param A number of maximum components
#' @param randomization if \code{TRUE} computes p-values. The harmonic mean
#' is considered across folds and times.
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return returns the mean performance across all folds and all repeats for
#' each components
#' @export
#' @importFrom caret createFolds
#' @importFrom psych harmonic.mean


crossValidation <- function(X, Y, K = 10, method = "cv", times = 1,seed = 123,
                   A = 5, randomization = FALSE, ...){


  if(method == "cv"){times = 1}
  set.seed(seed)

  res_out <- array(NA, dim = c(times, K, A))
  pv <- array(NA, dim = c(times, K, A))
  pv_adj <- array(NA, dim = c(times, K, A))

  tt <- 1
  while(times >= tt){

    flds <-  createFolds(c(1:dim(X)[1]), k = K)
    for(k in seq(K)){
      Xk <- X[-flds[[k]],]
      Yk <- Y[-flds[[k]]]

      a<-1
      while(a <= A){
          out_mcc <- mccTest(X = Xk, Y = Yk, A = a, randomization = randomization,...)

          res_out[tt, k, a] <- out_mcc$test

          if(randomization){
            pv[tt, k, a] <- out_mcc$pv
            pv_adj[tt, k, a] <- out_mcc$pv_adj
          }

          a <- a + 1
        }

      #end k folds
      }

      #end times
    tt <- tt+1
    }

  if(randomization){
 res <- list(pv =  apply(pv, 3, function(x) harmonic.mean(x))[1,],
             pv_adj = apply(pv_adj, 3, function(x) harmonic.mean(x))[1,],
             res_out = apply(res_out, 3, mean))
  }else{
    res <- apply(res_out, 3, mean)
  }
    return(res)
  }



