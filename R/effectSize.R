#' @title effect size PLSc
#' @description compute effect size based on MCC or predictive scores
#' @usage effectSize(X, Y, A, measure, randomization, ...)
#' @param X data matrix where columns represent the \eqn{p} classes and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and
#' rows the \eqn{n} observations.
#' @param A number of components
#' @param measure now \code{"MCC"} and \code{"score"} are implemented.
#' @param randomization compute p-values?
#' @param ... Futher parameters.
#' @author Angela Andreella
#' @return Returns the corresponding pvalues
#' @importFrom compositions clrInv
#' @importFrom compositions clr
#' @importFrom stats sd
#' @export

effectSize <- function(X, Y, A, measure = "score", randomization = FALSE, ...){

  if(measure == "mcc"){

    out <- mccTest(X = X, Y = Y, A = A, randomization = randomization,...)


  }

  if(measure == "score"){
    out <- scoreTest(X = X, Y = Y, A = A, randomization = randomization,...)
  }

  return(out)
}
