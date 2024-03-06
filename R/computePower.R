#' @title Power estimation
#' @description estimate power for a given sample size, alpha level and number
#' of score components.
#' @usage computePower(X, Y, A, n, seed = 123,
#' Nsim = 100, nperm = 200, alpha = 0.05,
#' test = "R2", Y.prob = FALSE, eps = 0.01, ...)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param A number of score components
#' @param n sample size
#' @param seed seed value
#' @param Nsim number of simulations
#' @param nperm number of permutations
#' @param alpha type I error
#' @param test type of test, one of \code{c("score", "mcc", "R2")}.
#' @param Y.prob Boolean value. Default @FALSE. IF @TRUE \code{Y} is a probability vector
#' @param eps Default 0.01. \code{eps} is used when \code{Y.prob = FALSE} to transform \code{Y} in a probability vector.
#' Default to "R2".
#' @param ... Futher parameters see \code{\link{PLSc}}
#' @author Angela Andreella
#' @return Returns the corresponding estimated power
#' @export
#' @importFrom progress progress_bar
#' @examples
#' \dontrun{
#' datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- computePower(X = datas$X, Y = datas$Y, A = 3, n = 20)
#' }



computePower <- function(X, Y, A, n, seed = 123,
                         Nsim = 100, nperm = 200, alpha = 0.05,
                         test = "R2", Y.prob = FALSE, eps = 0.01, ...){

  #Build the reference model PLS2c

  outPLS <- PLSc(X = X, Y = Y, A = A, Y.prob = Y.prob, eps = eps, ...)

  pw <- matrix(0, ncol = length(test), nrow = A)
  colnames(pw)<- test

  pb <- progress_bar$new(total = Nsim)
  for(i in seq(Nsim)){
    #sapply(seq(Nsim), function(x) {
    #Model the distribution of the X-data
    outsim <- sim_XY(out = outPLS, n = n, seed = 1234+i, A = A, ...)
    #Model the distribution of the Y-data
    Xsim <- outsim$X_H1
    Ysim <- outsim$Y_H1


    #Apply one test

    if(length(test) == 1){
   # if(test == "eigen"){
  #    pv <- eigenTest(X = Xsim, Y = Ysim, A = A, nperm = nperm, ...)
  #  }

    if(test == "mcc"){
      pv <- sapply(seq(A), function(x) mccTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                    randomization = TRUE, ...))
      pv <- data.frame(pv = unlist(as.matrix(pv)[1,]), pv_adjust = unlist(as.matrix(pv)[2,]))
    }
    if(test == "score"){
      pv <- sapply(seq(A), function(x) scoreTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                    randomization = TRUE,...))
      pv <- data.frame(pv = unlist(as.matrix(pv)[1,]), pv_adjust = unlist(as.matrix(pv)[2,]))
    }
      if(test == "R2"){
        pv <- sapply(seq(A), function(x) R2Test(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                                                   randomization = TRUE, ...))
        pv <- data.frame(pv = unlist(as.matrix(pv)[1,]), pv_adjust = unlist(as.matrix(pv)[2,]))
      }

    for(x in seq(A)){
      if(pv$pv_adj[x] <= alpha){pw[x] <- pw[x] + 1}
    }
    }else{

      #Apply more than one test.
      pv <- data.frame(NA)

   #   if("eigen" %in% test){
    #    pv <- cbind(pv, eigen = eigenTest(X = Xsim, Y = Ysim, A = A, nperm = nperm,
    #                    scaling = scaling, ...)[["pv_adj"]])
   #   }

      if("mcc" %in% test){
        pv_out <- sapply(seq(A), function(x) mccTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                                                     randomization = TRUE,...))

        pv_out <- data.frame(pv = unlist(as.matrix(pv_out)[1,]), pv_adjust = unlist(as.matrix(pv_out)[2,]))

        pv <- cbind(pv, mcc = pv_out$pv_adjust)
      }
      if("score" %in% test){
        pv_out <- sapply(seq(A), function(x) scoreTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                                                       randomization = TRUE,...))

        pv_out <- data.frame(pv = unlist(as.matrix(pv_out)[1,]), pv_adjust = unlist(as.matrix(pv_out)[2,]))
        pv <- cbind(pv, score = pv_out$pv_adjust)

      }
      if("R2" %in% test){
        pv_out <- sapply(seq(A), function(x) R2Test(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                                                    randomization = TRUE,...))

        pv_out <- data.frame(pv = unlist(as.matrix(pv_out)[1,]), pv_adjust = unlist(as.matrix(pv_out)[2,]))
        pv <- cbind(pv, score = pv_out$pv_adjust)

      }

      pv <- pv[,-1]

      for(x in seq(A)){
        for(y in seq(length(test))){
          if(pv[x,y] <= alpha){
            pw[x,y] <- pw[x,y] + 1
          }
        }

      }

    }



    pb$tick()
    Sys.sleep(1 / Nsim)
 #   })
}
  pw <- pw/Nsim

  return(pw)
}
