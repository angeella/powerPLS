#' @title sample size calculation
#' @description Compute optimal sample size calculation
#' @usage computeSampleSize(n, X, Y, A, post.transformation, alpha, beta,
#' nperm, Nsim, seed, scaling, test = "mcc",...)
#' @param n vector of sample sizes to consider
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} classes and
#' rows the \eqn{n} observations.
#' @param A number of latent components
#' @param post.transformation TRUE if you want to apply post transformation.
#' @param alpha level of type I error
#' @param beta level of type II error
#' @param Nsim number of simulations
#' @param nperm number of permutations
#' @param seed seed value
#' @param scaling type of scaling, one of
#' \code{c("auto-scaling", "pareto-scaling", "mean-centering")}
#' @param test type of test, one of \code{c("score", "mcc", "eigen")}.
#' Default to "mcc".
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns a dataframe that contains the estimated power for each
#' sample size and number of components considered
#' @export


computeSampleSize <- function(n, X, Y, A, post.transformation, alpha, beta,
                              nperm, Nsim, seed, scaling, test = "mcc",...){


  samplesize <-  sapply(seq(length(n)), function(x) computePower(X = X, Y = Y, A = A,
                                 post.transformation = post.transformation,
                                 n = n[x], nperm = nperm, Nsim = Nsim,
                                 scaling = scaling, alpha = alpha, test = test, ...))


  samplesize <- list(size = n,
                     power = samplesize,
                     A = A)

  out <- data.frame(Power = as.vector(t(samplesize$power)),
                    Size = rep(n ),
                    A =rep(seq(A), each =length(n)))


  return(out)
}
