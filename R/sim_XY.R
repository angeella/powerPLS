#' @title simulate X data matrix under the alternative hypothesis
#' @description simulate X data matrix under the alternative hypothesis.
#' @usage sim_XY(X, Y, A, nperm, A, scaling, post.transformation)
#' @param out output from PLS2c
#' @param n number of observations
#' @param A number of components
#' @param post.transformation TRUE if you want to perform post transformation.
#' @param seed fix seed
#' @author Angela Andreella
#' @return Returns a simulated X matrix under the alternative hypothesis.
#' @export
#' @importFrom simukde simulate_kde


sim_XY <- function(out, n, seed = 123, post.transformation, A){


  set.seed(seed)
#  out <- PLSc(X = X, Y = Y, A = A, scaling = scaling, post.transformation = post.transformation)
  if(post.transformation){
    M <- out$M
  }else{
    M <- A
  }

  T_score <- out$T.score
  X_loading <- out$X.loading
  Y_loading <- out$Y.loading
  B <- out$B

  if(A == M){
    T_scoreO <- T_score[,1:A]
  sim_TO <- sapply(seq(A), function(x) simulate_kde(x = T_scoreO[,x], n = n)$random.values)
    sim_TO <- sim_Tscore(Tscore = T_scoreO, n = n, seed = sample.int(1))

     T_sim <- sim_TO #T target

  }else{
    T_scoreO <- T_score[,1:M]
    if(is.null(ncol(T_score[,((M+1):A)]))){
      ncp <- 1
    }else{
      ncp <- ncol(T_score[,((M+1):A)])
    }

    if(is.null(ncol(T_score[,1:M]))){
      nco <- 1
    }else{
      nco <- ncol(T_score[,1:M])
    }
    T_scoreO <- matrix(T_score[,(1:M)], ncol = nco)
    T_scoreP <- matrix(T_score[,((M+1):A)], ncol = ncp)
    sim_TO <- unlist(sapply(seq(ncol(T_scoreO)), function(x) simulate_kde(x = T_scoreO[,x], n = n)$random.values))
    sim_TP <- unlist(sapply(seq(ncol(T_scoreP)), function(x) simulate_kde(x = T_scoreP[,x], n = n)$random.values))
    #sim_TO <- simulate_kde(x = T_scoreO, n = n)$random.values
   # sim_TP <-simulate_kde(x = T_scoreP, n = n)$random.values
    T_sim <- cbind(sim_TO, sim_TP) #T target
  }

  T_sim <- apply(T_sim, 2, function(x) scale(x, scale = FALSE))
  out1 <- svd(T_sim %*% t(T_score) %*% T_score)

  S <- (t(T_score) %*% T_score)
  T_new <- out1$u %*% t(out1$v) %*% diag(diag(S))^(1/2)

  E_pilot <- X - T_score %*% t(X_loading)

  E <- permuteIndex(E_pilot, times = n, by.row = TRUE, replace = TRUE)

  E_X <- (diag(dim(T_new)[1]) - T_new %*% solve(t(T_new) %*% T_new) %*% t(T_new)) %*% E

  X_H1 <- T_new %*% t(X_loading) + E_X


  U <- runif(1, 0, 1)
  Y_probH1 <- matrix(clrInv(X %*% B)$Y, ncol = 1)

  Y_H1 <- sapply(c(1:nrow(Y)), function(x) ifelse(Y_probH1[x] > runif(1, 0, 1), 0, 1))
 # Y_H1 <- apply(Y_H1, 1, which.max)

#  Y_H1 <- dummy(Y_H1)

  return(list(Y_H1 = Y_H1, X_H1 = X_H1))
}
