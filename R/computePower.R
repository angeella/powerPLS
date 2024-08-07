#' @title Power estimation
#' @description estimate power for a given sample size, alpha level and number
#' of score components.
#' @usage computePower(X, Y, A, n, seed = 123,
#' Nsim = 100, nperm = 200, alpha = 0.05,
#' scaling = "auto-scaling", test = "R2",
#' Y.prob = FALSE, eps = 0.01, post.transformation = TRUE,
#' fast=FALSE,transformation = "clr")
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
#' @param fast use fk_density from the FKSUM package for KDE. Default @FALSE.
#' @param transformation transformation used to map \code{Y} in probability data vector. The options are @ilr and @clr.
#' @author Angela Andreella
#' @return Returns the corresponding estimated power
#' @export
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
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
                         post.transformation = TRUE,
                         fast=FALSE,transformation = "clr") {

  if (any(!(test %in% c("R2", "mcc", "score")))) {
    stop("available tests are R2, mcc and score")
  }
  # Build the reference model PLS2c
  outPLS <- PLSc(X = X, Y = Y, A = A, Y.prob = Y.prob, eps = eps, scaling =
                 scaling, post.transformation = post.transformation,
                 transformation = "clr")

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    cl <- 2L # use 2 cores
  } else {
    cl <- parallel::detectCores() - 2 #not overload the computer
  }

  cl <- parallel::makeCluster(parallel::detectCores())

pw <- foreach(a = c(1:Nsim),.errorhandling = "remove")%dopar%{

  pw_sim <- matrix(0, ncol = length(test), nrow = A)

  outsim <- sim_XY(out = outPLS, n = n, seed = a, A = A,
                   post.transformation = post.transformation,fast=fast)

  Xsim <- outsim$X_H1
  rownames(outsim$Y_H1) <- NULL


  if(length(table(outsim$Y_H1))==1){
    Ysim <- as.numeric(outsim$Y_H1)
    if(Ysim[1]==0){
      Ysim[1] <- Ysim[1] +1
    }else{
      Ysim[1] <- Ysim[1] -1
    }

    }else{
    Ysim <- outsim$Y_H1[,2]
  }

  results <- list()

  if ("mcc" %in% test) {

    results$pv_mcc <- sapply(c(1:A), function(x){
                                mccTest(X = Xsim, Y = Ysim,
                                        nperm = nperm, A=x,
                                        randomization = TRUE,
                                        Y.prob = Y.prob, eps = eps,
                                        scaling = scaling,
                                        post.transformation = post.transformation)$pv_adj
                              })
  }
  if ("score" %in% test) {
    results$pv_score <- sapply(c(1:A), function(x){
                                  scoreTest(X = Xsim, Y = Ysim,
                                            nperm = nperm, A=x,
                                            randomization = TRUE,
                                            Y.prob = Y.prob, eps = eps,
                                            scaling = scaling,
                                            post.transformation = post.transformation)$pv_adj
                                })
  }
  if ("R2" %in% test) {
    results$pv_R2 <- sapply(c(1:A), function(x){
                               R2Test(X = Xsim, Y = Ysim,
                                      nperm = nperm, A=x,
                                      randomization = TRUE,
                                      Y.prob = Y.prob, eps = eps,
                                      scaling = scaling,
                                      post.transformation = post.transformation)$pv_adj
       })
  }

  pv_out <- data.frame(matrix(unlist(results), nrow = A))

  pw_sim<-ifelse(pv_out<=alpha, pw_sim +1, pw_sim)

  colnames(pw_sim) <- gsub("pv_", "", names(results))
  rownames(pw_sim) <- seq(A)
  pw_sim
  }
  parallel::stopCluster(cl)
  Nsim_final <- length(pw)
  if(Nsim != Nsim_final){
    warning(paste0("The power was calculated with "), Nsim_final, " simulations instead of ", Nsim)
  }
  pw <- Reduce('+', pw)/Nsim_final

  return(pw)
}
