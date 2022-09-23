
#Utils

#Fit Y
fitY <- function(X, B, Mm, s){

  Y.fitted = matrix(apply(compositions::clrInv(s*(X %*% B) + Mm), 1, function(x) which.max(x)),ncol = 1)
#  Y.fitted = matrix(Y.fitted, ncol = ncol(Y))
  Y.fitted <- ifelse(Y.fitted == 2, 1, 0)

  return(Y.fitted)
}

permuteIndex <- function(Y, by.row = TRUE, times, replace = FALSE){

  if(by.row){
    idx <- sample(seq(nrow(Y)), size = times, replace = replace)
    Y <- Y[idx,]
  }else{
    idx <- sample(seq(ncol(Y)), size =times, replace = replace)
    Y <- Y[,idx]
  }
return(Y)
}

similarityMatrix <- function(X, Y){

  out <- svd(X - Y)
  eigen_d <- sqrt(sum(log(out$d)^ 2))
  #Escoufier's R V Coefficient
  RV <- sum(diag(t(X) %*% Y %*% t(Y) %*% X))/sqrt(sum(diag(X %*% t(X)))^2 * sum(diag(X %*% t(X)))^2)

  return(data.frame(eigen_d = eigen_d, RV = RV))
}
