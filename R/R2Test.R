#' @title R2 test
#' @description Performs randomization test based on R2
#' @usage R2Test(X, Y, nperm = 100, A, randomization = FALSE, Y.prob = FALSE, eps = 0.01,
#' scaling = "auto-scaling", post.transformation = TRUE)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param nperm number of permutations. Default 100.
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
#' @export
#' @seealso The type of tests implemented: \code{\link{scoreTest}} \code{\link{mccTest}}.
#' @author Angela Andreella
#' @return List with the following objects: \code{pv}: raw p-value, \code{pv_adj}: adjusted p-value, \code{test} estimated statistical test.
#' @export
#' @references For the general framework of power analysis for PLS-based methods see:
#'
#' @examples
#' datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- R2Test(X = datas$X, Y = datas$Y, A = 1)
#' out


R2Test <- function(X, Y, nperm = 100, A, randomization = FALSE,
                   Y.prob = FALSE, eps = 0.01, scaling = "auto-scaling",
                   post.transformation = TRUE){


  out <- PLSc(X = X, Y = Y, A = A, scaling = scaling,
              post.transformation = post.transformation,eps = eps,
              Y.prob = Y.prob,transformation = "ilr")

  if(!Y.prob){

    if(is.null(dim(Y)) | ncol(as.matrix(Y))==1){
      Yf <- as.matrix(Y)

      if(!is.factor(Yf)){
        Yf <- as.factor(Yf)

      }
      levels(Yf) <- c(0,1)
      Yf <- model.matrix(~0+Yf)
    }

    #Transform to probability matrix
    Yf[which(Yf==0)]<-eps
    Yf[which(Yf==1)]<-1-(ncol(Yf)-1)*eps

    #Centered log ratio transform transformation
    P <- matrix(ilr(Yf), ncol = 1)
  }else{
    P <- Y
  }
  s <- sd(P)
  Yfitted = matrix(ilrInv(s*(X %*% out$B))[,3], ncol = 1)

  #observed R2
  r2_obs <- cor(Yfitted, P)[[1]]

  if(randomization){
    null_distr <- replicate(nperm-1, {

      idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
      Xkp <- X[idx,]

      out <- PLSc(X = Xkp, Y = Y, A = A, scaling = scaling,
                  post.transformation = post.transformation,eps = eps,
                  Y.prob = Y.prob,transformation = "ilr")


      Yfitted =matrix(ilrInv(s*(Xkp %*% out$B))[,3], ncol = 1)


      cor(Yfitted, P)[[1]]

    })
    null_distr <-c(r2_obs, null_distr)
    pv <- mean(null_distr >= r2_obs)

  }else{
    pv <- NA
  }



  return(list(pv = pv,
              pv_adj = min(pv*A,1), test = r2_obs))

}
