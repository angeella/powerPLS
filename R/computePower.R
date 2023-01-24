#' @title power computations
#' @description Performs randomization test based on eigenvalues
#' @usage computePower(X, Y, A, nperm, A, scaling, post.transformation)
#' @param X data matrix where columns represent the \eqn{p} variables and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} classes and rows the \eqn{n} observations.
#' @param nperm number of permutations
#' @param A number of components
#' @param post.transformation TRUE if you want to apply post transformation.
#' @param scaling type of scaling.
#' @param n sample size
#' @param seed seed value
#' @param Nsim number of simulations
#' @param nperm number of permutations
#' @param alpha type I error
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns the corresponding pvalues
#' @export
#' @importFrom progress progress_bar




computePower <- function(X, Y, A,
                          post.transformation, scaling, n, seed = 123,
                          Nsim, nperm, alpha, eps = 0.01, ...){

  #Build the reference model PLS2c

  outPLS <- PLSc(X = X, Y = Y, A = A,
              scaling = scaling, eps = eps,
              post.transformation = post.transformation, ...)

  pw <- 0
  pb <- progress_bar$new(total = Nsim)
  for(i in seq(Nsim)){
    #Model the distribution of the X-data
    outsim <- sim_XY(out = outPLS, n = n, seed = seed,
                     post.transformation = post.transformation, A = A)
    #Model the distribution of the Y-data
    Xsim <- outsim$X_H1
    Ysim <- outsim$Y_H1

    pv <- eigenTest(X = Xsim, Y = Ysim, A = A, nperm = nperm,
                    scaling = scaling, ...)



    if(pv$pv[A] <= alpha){pw <- pw + 1}
    pb$tick()
    Sys.sleep(1 / Nsim)
  }

  pw <- sum(pw)/Nsim

  return(pw)
}
