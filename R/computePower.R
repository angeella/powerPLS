#' @title power computations
#' @description Performs randomization test based on eigenvalues
#' @usage computePower(X, Y, A,
#' post.transformation = TRUE, scaling, n, seed = 123,
#' Nsim, nperm, alpha, eps = 0.01, test = "eigen", ...)
#' @param X data matrix where columns represent the \eqn{p} variables and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} classes and rows the \eqn{n} observations.
#' @param A number of components
#' @param post.transformation TRUE if you want to apply post transformation.
#' @param scaling type of scaling.
#' @param n number of observations, vector of length two (first values is the
#' number of observations for the first classe, the second for the second class)
#' @param seed seed value
#' @param Nsim number of simulations
#' @param nperm number of permutations
#' @param alpha type I error
#' @param eps see details
#' @param test type of test
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns the corresponding pvalues
#' @export
#' @importFrom progress progress_bar




computePower <- function(X, Y, A,
                          post.transformation = TRUE, scaling, n, seed = 123,
                          Nsim, nperm, alpha, eps = 0.01, test = "eigen", ...){

  #Build the reference model PLS2c

  outPLS <- PLSc(X = X, Y = Y, A = A,
              scaling = scaling, eps = eps,
              post.transformation = post.transformation, ...)

  pw <- rep(0, A)
  pb <- progress_bar$new(total = Nsim)
  for(i in seq(Nsim)){
    #sapply(seq(Nsim), function(x) {
    #Model the distribution of the X-data
    outsim <- sim_XY(out = outPLS, n = n, seed = 1234+i,
                     post.transformation = post.transformation, A = A)
    #Model the distribution of the Y-data
    Xsim <- outsim$X_H1
    Ysim <- outsim$Y_H1

    if(test == "eigen"){
      pv <- eigenTest(X = Xsim, Y = Ysim, A = A, nperm = nperm,
                      scaling = scaling, ...)
    }

    if(test == "mcc"){
      pv <- sapply(seq(A), function(x) mccTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                    scaling = scaling, randomization = TRUE, eps = eps,
                    post.transformation = post.transformation,...))
      pv <- data.frame(pv = unlist(as.matrix(pv)[1,]), pv_adjust = unlist(as.matrix(pv)[2,]))
    }
    if(test == "score"){
      pv <- sapply(seq(A), function(x) scoreTest(X = Xsim, Y = Ysim[,2], A = x, nperm = nperm,
                    scaling = scaling, randomization = TRUE, eps = eps,
                    post.transformation = post.transformation,...))
      pv <- data.frame(pv = unlist(as.matrix(pv)[1,]), pv_adjust = unlist(as.matrix(pv)[2,]))
    }


    for(x in seq(A)){
      if(pv$pv_adj[x] <= alpha){pw[x] <- pw[x] + 1}
    }

    pb$tick()
    Sys.sleep(1 / Nsim)
 #   })
}
  pw <- pw/Nsim

  return(pw)
}
