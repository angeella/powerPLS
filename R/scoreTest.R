#' @title Score test
#' @description Performs randomization test based on predictive score vector
#' @usage scoreTest(X, Y, nperm = 200, A, randomization = FALSE,
#' Y.prob = FALSE, eps = 0.01, scaling = "auto-scaling",
#' post.transformation = TRUE)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param nperm number of permutations. Default 200.
#' @param A number of score components
#' @param randomization Boolean value. Default @FALSE. If @TRUE the permutation p-value is computed
#' @param Y.prob Boolean value. Default @FALSE. IF @TRUE \code{Y} is a probability vector
#' @param eps Default 0.01. \code{eps} is used when \code{Y.prob = FALSE} to transform \code{Y} in a probability vector
#' @param scaling type of scaling, one of
#' \code{c("auto-scaling", "pareto-scaling", "mean-centering")}. Default @auto-scaling
#' @param post.transformation Boolean value. @TRUE if you want to apply post transformation. Default @TRUE
#' @author Angela Andreella
#' @return Returns a list with the corresponding statistical tests,
#' raw and adjusted p-values
#' @importFrom compositions ilr
#' @importFrom stats cor
#' @importFrom stats t.test
#' @importFrom stats var
#' @export
#' @seealso The type of tests implemented: \code{\link{mccTest}} \code{\link{R2Test}}.
#' @author Angela Andreella
#' @return List with the following objects: \code{pv}: raw p-value, \code{pv_adj}: adjusted p-value, \code{test} estimated statistical test.
#' @export
#' @references For the general framework of power analysis for PLS-based methods see:
#'
#' @examples
#' datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- scoreTest(X = datas$X, Y = datas$Y, A = 1)
#' out





scoreTest <- function(X, Y, nperm = 200, A, randomization = FALSE,
                      Y.prob = FALSE, eps = 0.01, scaling = "auto-scaling",
                      post.transformation = TRUE){

  out <- PLSc(X = X, Y = Y, A = A, scaling = scaling,
              post.transformation = post.transformation,
              eps = eps, Y.prob = Y.prob, transformation = "clr")


  T_score <- out$T_score

  if(!is.na(out$M)){
    M <- out$M
  }else{
    M <- A
  }

  if(A!=M){
    Tp <- T_score[,(M+1):A]
  }else{
    Tp <- T_score
  }

  lev <- unique(as.vector(Y))

  if(length(lev)!=2){stop("Y must be a binary variable")}

  if((sum(Y==lev[1])==1)|(sum(Y==lev[2])==1)){

  if(is.null(dim(Tp))){
      effect_obs <- t.test(Tp[Y == lev[1]]- Tp[Y == lev[2]], var.equal = FALSE)$statistic
  }else{
      effect_obs <- t.test(Tp[Y == lev[1],]- Tp[Y == lev[2],], var.equal = FALSE)$statistic
  }

  }else{

  if(is.null(dim(Tp))){
    effect_obs <- t.test(Tp[Y == lev[1]], Tp[Y == lev[2]], var.equal = FALSE)$statistic
  }else{
    effect_obs <- t.test(Tp[Y == lev[1],], Tp[Y == lev[2],], var.equal = FALSE)$statistic
  }
  }
  if(randomization){

    null_distr <- replicate(nperm-1,{

    idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
    Xkp <- X[idx,]

    out <- PLSc(X = Xkp, Y = Y, A = A, scaling = scaling,
                post.transformation = post.transformation,
                eps = eps, Y.prob = Y.prob, transformation = "clr")


      T_score <- out$T_score
      T_score

      if(!is.na(out$M)){
        M <- out$M
      }else{
        M <- A
      }

      if(A!=M){
        Tp <- T_score[,(M+1):A]
      }else{
        Tp <- T_score
      }

      if((sum(Y==lev[1])==1)|(sum(Y==lev[2])==1)){

        if(is.null(dim(Tp))){
          effect_obs <- t.test(Tp[Y == lev[1]]- Tp[Y == lev[2]], var.equal = FALSE)$statistic
        }else{
          effect_obs <- t.test(Tp[Y == lev[1],]- Tp[Y == lev[2],], var.equal = FALSE)$statistic
        }

      }else{

        if(is.null(dim(Tp))){
          effect_obs <- t.test(Tp[Y == lev[1]], Tp[Y == lev[2]], var.equal = FALSE)$statistic
        }else{
          effect_obs <- t.test(Tp[Y == lev[1],], Tp[Y == lev[2],], var.equal = FALSE)$statistic
        }
      }



    })

    null_distr <-abs(c(effect_obs, null_distr))
    pv <- mean(null_distr >= effect_obs)
  }else{
    pv <- NA
  }


  return(list(pv = pv,
              pv_adj = min(pv*A,1), test = effect_obs))
}
