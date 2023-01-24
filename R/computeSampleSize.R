#' @title sample size calculation
#' @description Compute optimal sample size calculation
#' @usage computeSampleSize(n, A, nperm, A, scaling, post.transformation)
#' @param n two dimensional vector giving the minimum and maximum sample size
#' to be considered
#' @param X data matrix where columns represent the \eqn{p} variables and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} classes and rows the \eqn{n} observations.
#' @param A number of components
#' @param post.transformation TRUE if you want to apply post transformation.
#' @param alpha level of type I error
#' @param beta level of type II error
#' @param nperm number of permutations
#' @param Nsim number of simulations
#' @param seed seed value
#' @param scaling type of scaling.
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns the corresponding pvalues


#add power?

computeSampleSize <- function(n, X, Y, A, post.transformation, alpha, beta,
                              nperm, Nsim, seed, scaling, ...){


  for(i in seq(min(n), max(n))){

  samplesize[i] <-  computePower(X = X, Y = Y, A = A,
                                 post.transformation = post.transformation,
                                 n = i, nperm = nperm, Nsim = Nsim, scaling = scaling, ...)
  }


  n <- which(seq(min(n), max(n))[samplesize == 1-beta])

  return(n)
}
