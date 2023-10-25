#' @title score test
#' @description Performs randomization test based on scores
#' @usage scoreTest(X, Y, nperm, A, randomization, ...)
#' @param X data matrix where columns represent the \eqn{p} classes and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and
#' rows the \eqn{n} observations.
#' @param nperm number of permutations
#' @param A number of components
#' @param randomization compute p-values?
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns a list with the corresponding statistical tests,
#' raw and adjusted p-values
#' @importFrom compositions clrInv
#' @importFrom compositions clr
#' @importFrom stats sd
#' @importFrom stats var
#' @export




scoreTest <- function(X, Y, nperm = 100, A, randomization = FALSE, scaling,...){

#  out <- PLSc(X = X, Y = Y, A = A, ...)
  out <- PLSc(X = X, Y = Y, A = A, scaling = scaling,
              post.transformation = TRUE,eps = 0.01,Y.prob = FALSE)


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

  if(is.null(dim(Tp))){
    effect_obs <- abs(mean(Tp[Y == lev[1]]) -  mean(Tp[Y == lev[2]]))/
      sqrt((var(Tp[Y == lev[1]]) +  var(Tp[Y == lev[2]]))/2)
  }else{
    effect_obs <- abs(mean(Tp[Y == lev[1],]) -  mean(Tp[Y == lev[2],]))/
      sqrt((var(as.vector(Tp[Y == lev[1],])) +  var(as.vector(Tp[Y == lev[2],])))/2)
  }

  if(randomization){
    pv <- 0
    for(j in seq(nperm-1)){
      idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
      Xkp <- X[idx,]

      #out <- PLSc(X = Xkp, Y = Y, A = A, ...)
     out <- PLSc(X = Xkp, Y = Y, A = A, scaling = scaling,post.transformation = TRUE,eps = 0.01,Y.prob = FALSE)
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
      lev <- unique(Y)
      if(is.null(dim(Tp))){
        effect_p <- abs(mean(Tp[Y == lev[1]]) -  mean(Tp[Y == lev[2]]))/
          sqrt((var(Tp[Y == lev[1]]) +  var(Tp[Y == lev[2]]))/2)
      }else{
        effect_p <- abs(mean(Tp[Y == lev[1],]) -  mean(Tp[Y == lev[2],]))/
          sqrt((var(as.vector(Tp[Y == lev[1],])) +  var(as.vector(Tp[Y == lev[2],])))/2)
      }

      if(is.na(effect_p)){effect_p <- effect_obs}
      if(effect_p >= effect_obs){pv <- pv+1}

    }

    pv <- (pv+1)/nperm
  }else{
    pv <- NA
  }


  return(list(pv = pv,
              pv_adj = min(pv*A,1), test = effect_obs))
}
