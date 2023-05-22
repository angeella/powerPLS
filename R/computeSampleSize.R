#' @title sample size calculation
#' @description Compute optimal sample size calculation
#' @usage computeSampleSize(n, X, Y, A, post.transformation, alpha, beta,
#' nperm, Nsim, seed, scaling, test = "eigen",...)
#' @param n list, where each objects is a vector of two dimensions representing the
#' number of observations for each class.
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
#' @param test type of test
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns the corresponding pvalues
#' @export


computeSampleSize <- function(n, X, Y, A, post.transformation, alpha, beta,
                              nperm, Nsim, seed, scaling, test = "eigen",...){


  samplesize <-  sapply(seq(length(n)), function(x) computePower(X = X, Y = Y, A = A,
                                 post.transformation = post.transformation,
                                 n = n[x][[1]], nperm = nperm, Nsim = Nsim,
                                 scaling = scaling, alpha = alpha, test = test, ...))


  samplesize <- list(size = n,
                     power = samplesize,
                     A = A)
  classprop <- sapply(seq(length(n)), function(x) n[[x]][1]/n[[x]][2])
  n <- sapply(seq(length(n)), function(x) sum(n[[x]]))
  classprop <-
  out <- data.frame(Power = as.vector(t(samplesize$power)),
                    Size = rep(n),
                    A =rep(seq(A), each =length(n)),
                    classprop = rep(classprop))

  return(out)
}
