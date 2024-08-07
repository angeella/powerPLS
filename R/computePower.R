#' @title Power estimation
#' @description estimate power for a given sample size, alpha level and number
#' of score components.
#' @usage computePower(X, Y, A, n, seed = 123,
#' Nsim = 100, nperm = 200, alpha = 0.05,
#' scaling = "auto-scaling", test = "R2",
#'Y.prob = FALSE, eps = 0.01, post.transformation = TRUE,...)
#' @param X data matrix where columns represent the \eqn{p} variables and
#' rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the two classes and
#' rows the \eqn{n} observations.
#' @param A number of score components
#' @param n sample size
#' @param seed seed value
#' @param Nsim number of simulations
#' @param nperm number of permutations
#' @param alpha type I error
#' @param scaling type of scaling, one of
#' \code{c("auto-scaling", "pareto-scaling", "mean-centering")}. Default @auto-scaling
#' @param test type of test, one of \code{c("score", "mcc", "R2")}.
#' @param Y.prob Boolean value. Default @FALSE. IF @TRUE \code{Y} is a probability vector
#' @param eps Default 0.01. \code{eps} is used when \code{Y.prob = FALSE} to transform \code{Y} in a probability vector.
#' Default to "R2".
#' @param post.transformation Boolean value. @TRUE if you want to apply post transformation. Default @TRUE
#' @param ... Futher parameters see \code{\link{PLSc}}
#' @author Angela Andreella
#' @return Returns the corresponding estimated power
#' @export
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#' @examples
#' \dontrun{
#' datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
#' out <- computePower(X = datas$X, Y = datas$Y, A = 3, n = 20, test = "R2")
#' }
#' @references
#'
#' Andreella, A., Finos, L., Scarpa, B. and Stocchero, M. "Towards a power analysis for PLS-based methods" 	arXiv:2403.10289 stat.ME.
#'



computePower <- function(X, Y, A, n, seed = 123,
                         Nsim = 100, nperm = 200, alpha = 0.05,
                         scaling = "auto-scaling",
                         test = "R2", Y.prob = FALSE, eps = 0.01,
                         post.transformation = TRUE,...) {

  if (any(!(test %in% c("R2", "mcc", "score")))) {
    stop("available tests are R2, mcc and score")
  }

  # Build the reference model PLS2c
  outPLS <- PLSc(X = X, Y = Y, A = A, Y.prob = Y.prob, eps = eps, scaling = scaling, post.transformation = post.transformation,...)

  # Funzione per eseguire una singola simulazione
  simulate_once <- function(i,...) {
    # Model the distribution of the pilot data
    outsim <- sim_XY(out = outPLS, n = n, seed = 1234 + i, A = A, post.transformation = post.transformation)

    Xsim <- outsim$X_H1
    Ysim <- outsim$Y_H1

    results <- list()

    if ("mcc" %in% test) {
      results$pv_mcc <- foreach(x = seq(A), .combine = 'cbind', .packages = c("powerPLS")) %dopar% {
        mccTest(X = Xsim, Y = Ysim[, 2], A = x, nperm = nperm, randomization = TRUE, Y.prob = Y.prob,eps = eps)
      }
    }
    if ("score" %in% test) {
      results$pv_score <- foreach(x = seq(A), .combine = 'cbind', .packages = c("powerPLS")) %dopar% {
        scoreTest(X = Xsim, Y = Ysim[, 2], A = x, nperm = nperm, randomization = TRUE, Y.prob = Y.prob, eps = eps)
      }
    }
    if ("R2" %in% test) {
      results$pv_R2 <- foreach(x = seq(A), .combine = 'cbind', .packages = c("powerPLS")) %dopar% {
        R2Test(X = Xsim, Y = Ysim[, 2], A = x, nperm = nperm, randomization = TRUE, Y.prob = Y.prob, eps = eps)
      }
    }

    pv_out <- data.frame(matrix(unlist(results), ncol = length(test), nrow=A))

    colnames(pv_out) <- test
    rownames(pv_out) <- seq(A)

    pw_sim <- matrix(0, ncol = length(test), nrow = A)
    colnames(pw_sim) <- test

    for (x in seq(A)) {
      for (y in seq(length(test))) {
        if (pv_out[x, y] <= alpha) {
          pw_sim[x, y] <- pw_sim[x, y] + 1
        }
      }
    }

    return(pw_sim)
  }

  # Numero di core da utilizzare

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    num_cores <- 2L
  } else {
    # use all cores in devtools::test()
    num_cores <- parallel::detectCores() -2
  }

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Esegui le simulazioni in parallelo usando foreach
  pw <- foreach(i = seq(Nsim), .combine = "+", .packages = c("powerPLS", "foreach")) %dopar% {
    simulate_once(i)
  }
  pw/Nsim

  stopCluster(cl)


  return(pw)
}
