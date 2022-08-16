#' @title sample size calculation
#' @description Compute optimal sample size calculation
#' @usage computeSampleSize(d, nperm, A, Nsim, post.transformation, scaling, seed, nperm, alpha = 0.05, eps = 0.01)
#' @param d effect size value or list containing X and Y data
#' @param nperm number of permutations
#' @param A number of components
#' @param Nsim number of simulations
#' @param power power desired
#' @param N initial number of observations
#' @param post.transformation TRUE if you want to apply post transformation.
#' @param scaling type of scaling
#' @param seed seed
#' @param alpha level of significance
#' @param eps epsilon value
#' @export
#' @author Angela Andreella
#' @return Returns the corresponding pvalues


#add power?

computeSampleSize <- function(d, nperm, A, Nsim, power, N = NULL, post.transformation, scaling, seed, alpha = 0.05, eps = 0.01){


  if(is.list(d)){
    #Retrospective power analysis
    X <- d$X
    Y <- d$Y
  }else{
    #Prospective power analysis

    #TODO: simulate data given d.
  }

  cond <- TRUE
  if(is.null(N)){N <- nrow(X)}
  while(cond){

    id <- sample(c(1:nrow(X)), N, replace = TRUE)

    Xsub <- X[id,]
    Ysub <- Y[id,]

    out <- computePower(X = Xsub, Y = Ysub, A = A, scaling = scaling,
                        post.transformation = post.transformation,
                        seed, n = nrow(Xsub),
                        Nsim = Nsim, nperm = nperm, alpha = alpha, eps = eps)

    if(out >= power){cond <- FALSE}else{N <- N + 1}
  }




return(out)
}
