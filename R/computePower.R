#' @title power computations
#' @description Performs randomization test based on eigenvalues
#' @usage computePower(X, Y, A, nperm, A, scaling, post.transformation)
#' @param X data matrix where columns represent the \eqn{p} classes and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and rows the \eqn{n} observations.
#' @param A number of components
#' @param X.scaling scaling X?
#' @param post.transformation TRUE if you want to apply post transformation.
#' @param seed fix seed
#' @param n sample size
#' @param Nsim number of simulations
#' @param nperm number of permutations
#' @param alpha level of significance
#' @param eps epsilon value used to transformed Y in probability data matrix.
#' @author Angela Andreella
#' @return Returns the corresponding pvalues
#' @export
#' @importFrom progress progress_bar




computePower <- function(X, Y, A, scaling,
                          post.transformation, seed, n = nrow(X),
                          Nsim = 100, nperm = 100, alpha = 0.05, eps = 0.01){

  #Build the reference model PLS2c

  outPLS <- PLSc(X = X, Y = Y, A = A,
              scaling = scaling, eps = eps,
              post.transformation = post.transformation)

  pw <- 0
  pb <- progress_bar$new(total = Nsim)
  for(i in seq(Nsim)){
    #Model the distribution of the X Y-data
    outsim <- sim_XY(outPLS, n = n, seed = sample.int(1),
                     post.transformation = post.transformation, A = A)

    Xsim <- outsim$X_H1
    Ysim <- outsim$Y_H1

    Ysim <- cbind(Ysim, ifelse(Ysim == 1, 0, 1))

    pv <- eigenTest(X = Xsim, Y = Ysim, A = A, nperm = nperm)



    if(pv$pv_adj[A] <= alpha){pw <- pw + 1}
    pb$tick()
    Sys.sleep(1 / Nsim)
  }

  pw <- sum(pw)/Nsim

  return(pw)
}
