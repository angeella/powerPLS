#' @title dQ2 test
#' @description Performs permutation-based test based on dQ2
#' @usage dQ2Test(X, Y, nperm = 200, A, randomization = FALSE,
#' Y.prob = FALSE, eps = 0.01, scaling = 'auto-scaling',
#' post.transformation = TRUE, class = 1, cross.validation = FALSE, ...)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param nperm number of permutations. Default to 200.
#' @param A number of score components
#' @param randomization Boolean value. Default to \code{FALSE}. If \code{TRUE} the permutation p-value is computed
#' @param Y.prob Boolean value. Default \code{FALSE}. IF \code{TRUE} \code{Y} is a probability vector
#' @param eps Default 0.01. \code{eps} is used when \code{Y.prob = FALSE} to transform \code{Y} in a probability vector
#' @param scaling Type of scaling, one of
#' \code{c('auto-scaling', 'pareto-scaling', 'mean-centering')}. Default 'auto-scaling'.
#' @param post.transformation Boolean value. \code{TRUE} if you want to apply post transformation. Default \code{TRUE}
#' @param class Numeric value. Specifiy the reference class. Default \code{1}
#' @param cross.validation Boolean value. Default \code{FALSE}. \code{TRUE} if you want to compute the observed test statistic by Nested cross-validation
#' @param ... additional arguments related to \code{cross.validation}. See \code{\link{repeatedCV_test}}
#' @author Angela Andreella
#' @importFrom compositions ilr
#' @importFrom stats cor
#' @export
#' @seealso Other test statistics implemented: \code{\link{mccTest}}, \code{\link{scoreTest}},
#' \code{\link{sensitivityTest}}, \code{\link{specificityTest}},\code{\link{AUCTest}}, \code{\link{R2Test}},
#' \code{\link{FMTest}}, \code{\link{F1Test}}.
#' @author Angela Andreella
#' @return List with the following objects:
#' \describe{
#'   \item{pv}{raw p-value. It equals \code{NA} if \code{randomization = FALSE}}
#'   \item{pv_adj}{adjusted p-value. It equals \code{NA} if \code{randomization = FALSE}}
#'   \item{test}{estimated test statistic}
#' }
#' @export
#' @references For the general framework of power analysis for PLS-based methods see:
#'
#' Andreella, A., Fino, L., Scarpa, B., & Stocchero, M. (2024). Towards a power analysis for PLS-based methods. arXiv preprint \url{https://arxiv.org/abs/2403.10289}.
#' @examples
#' datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 1)
#' out <- dQ2Test(X = datas$X, Y = datas$Y, A = 1)
#' out




dQ2Test <- function(X, Y, nperm = 200, A, randomization = FALSE, Y.prob = FALSE, eps = 0.01,
                    scaling = "auto-scaling", post.transformation = TRUE, class = 1, cross.validation = FALSE,
                    ...) {

  out <- PLSc(X = X, Y = Y, A = A, transformation = "clr", scaling = scaling, post.transformation = post.transformation,
              eps = eps, Y.prob = Y.prob)

  # Fitted from pilot data
  rownames(out$Y_fitted) <- NULL
  if (is.null(dim(out$Y_fitted))) {
    Y_fitted <- as.factor(out$Y_fitted)

  } else {
    Y_fitted <- as.factor(out$Y_fitted[, 2])

  }

  # Observed one
  Yf <- as.factor(Y)

  # check levels
  levels(Yf) <- c(0, 1)
  levels(Y_fitted) <- c(0, 1)

  # dQ2 observed
  if (cross.validation) {
    dQ2_obs <- repeatedCV_test(data = X, labels = Y, A = A, test_type = "dQ2Test", ...)
  } else {
    dQ2_obs <- dQ2(Yf, Y_fitted, class = class)
  }


  if (randomization) {
    null_distr <- replicate(nperm - 1, {
      # Permute rows
      idx <- sample(seq(nrow(X)), nrow(X), replace = FALSE)
      Xkp <- X[idx, ]

      # Compute Y fitted
      out <- PLSc(X = Xkp, Y = Y, A = A, transformation = "clr", scaling = scaling,
                  post.transformation = post.transformation, eps = eps, Y.prob = Y.prob)

      rownames(out$Y_fitted) <- NULL

      # Compute permuted dQ2
      Y_fitted <- as.factor(out$Y_fitted[, 2])
      levels(Y_fitted) <- c(0, 1)
      dQ2(Yf, Y_fitted, class = class)



    })

    null_distr <- c(dQ2_obs, null_distr)
    pv <- mean(null_distr >= dQ2_obs)


  } else {
    pv <- NA
  }


  return(list(pv = pv, pv_adj = min(pv * A, 1), test = dQ2_obs))
}
