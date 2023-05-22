#' @title simulate data matrix under the alternative hypothesis
#' @description simulate data matrix under the alternative hypothesis.
#' @usage sim_XY(out, n, seed = 123, post.transformation, A)
#' @param out output from PLS2c
#' @param n number of observations, vector of length two (first values is the
#' number of observations for the first classe, the second for the second class)
#' @param seed fix seed
#' @param post.transformation TRUE if you want to apply post transformation.
#' @param A number of components
#' @author Angela Andreella
#' @return Returns a simulated matrix under the alternative hypothesis.
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

  T_score <- out$T_score
  X_loading <- out$X_loading
  Y_loading <- out$Y_loading
  B <- out$B
  X <- out$X
  Y <- out$Y
  idx_class1 <- which(apply(as.matrix(Y), 1, which.max)==1)
  idx_class2 <- which(apply(as.matrix(Y), 1, which.max)==2)
  n1 <- n[1]
  n2 <- n[2]
  n <- sum(n)
  #First class

  T_score1 <- T_score[idx_class1,]
  if(A == M){
    T_scoreO <- T_score1
    sim_TO <- sapply(seq(A), function(x) simulate_kde(x = T_scoreO[,x], n = n1)$random.values)
    T_sim <- sim_TO #T target

  }else{
    T_scoreO <- T_score1[,1:M]
    if(is.null(ncol(T_score1[,((M+1):A)]))){
      ncp <- 1
    }else{
      ncp <- ncol(T_score1[,((M+1):A)])
    }

    if(is.null(ncol(T_score1[,1:M]))){
      nco <- 1
    }else{
      nco <- ncol(T_score1[,1:M])
    }
    T_scoreO <- matrix(T_score1[,(1:M)], ncol = nco)
    T_scoreP <- matrix(T_score1[,((M+1):A)], ncol = ncp)
    sim_TP <- unlist(sapply(seq(ncol(T_scoreP)), function(x) simulate_kde(x = T_scoreP[,x], n = n1)$random.values))
    sim_TO <- unlist(sapply(seq(ncol(T_scoreO)), function(x) simulate_kde(x = T_scoreO[,x], n = n1)$random.values))
    T_sim <- cbind(sim_TO, sim_TP) #T target

  }

  T_sim1 <- T_sim
  #Second class

  T_score1 <- T_score[idx_class2,]
  if(A == M){
    T_scoreO <- T_score1
    sim_TO <- sapply(seq(A), function(x) simulate_kde(x = T_scoreO[,x], n = n2)$random.values)
    T_sim <- sim_TO #T target

  }else{
    T_scoreO <- T_score1[,1:M]
    if(is.null(ncol(T_score1[,((M+1):A)]))){
      ncp <- 1
    }else{
      ncp <- ncol(T_score1[,((M+1):A)])
    }

    if(is.null(ncol(T_score1[,1:M]))){
      nco <- 1
    }else{
      nco <- ncol(T_score1[,1:M])
    }
    T_scoreO <- matrix(T_score1[,(1:M)], ncol = nco)
    T_scoreP <- matrix(T_score1[,((M+1):A)], ncol = ncp)
    sim_TP <- unlist(sapply(seq(ncol(T_scoreP)), function(x) simulate_kde(x = T_scoreP[,x], n = n2)$random.values))
    sim_TO <- unlist(sapply(seq(ncol(T_scoreO)), function(x) simulate_kde(x = T_scoreO[,x], n = n2)$random.values))
    T_sim <- cbind(sim_TO, sim_TP) #T target

  }

  T_sim <- rbind(T_sim1, T_sim)

  out1 <- svd(T_sim %*% t(T_score) %*% T_score)

  S <- (t(T_score) %*% T_score)
  if(length(S)==1){
    T_new <- out1$u %*% t(out1$v) %*% diag(S)^(1/2)
  }else{
    T_new <- out1$u %*% t(out1$v) %*% diag(diag(S))^(1/2)
  }


  E_pilot <- X - T_score %*% t(X_loading)

  E <- permuteIndex(E_pilot, times = n, by.row = TRUE, replace = TRUE)
  E <- as.matrix(E)

  E_X <- (diag(dim(T_new)[1]) - T_new %*% solve(t(T_new) %*% T_new) %*% t(T_new)) %*% E

  X_H1 <- T_new %*% t(X_loading) + E_X


  Y_H1 <- fitY(X = X_H1, B = B, Mm = 0, s = 1)
 # Y_H1 <- Y
#  eps <- min(unique(Y_H1))
#  Y_H1[which(Y_H1==eps)]<-0
#  Y_H1[which(Y_H1==1-eps)]<-1

  return(list(Y_H1 = Y_H1, X_H1 = X_H1))
}
