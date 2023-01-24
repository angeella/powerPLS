
#Utils

#Fit Y
fitY <- function(X, B, Mm, s){

  Y.fitted = matrix(apply(compositions::clrInv(s*(X %*% B) + Mm), 1, function(x) which.max(x)),ncol = 1)
  Y.fitted = as.factor(Y.fitted)
  Y.fitted = model.matrix(~0+Y.fitted)

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

  tr <- function(X){

    return(sum(diag(X)))
  }
  n <- nrow(X)
  XX <- crossprod(X)
  YY <- crossprod(Y)

  SXY <- (1/(n-1)) * t(XX) %*% YY
  SYX <- (1/(n-1)) * t(YY) %*% XX
  SXX <- (1/(n-1)) * t(XX) %*% XX
  SYY <- (1/(n-1)) * t(YY) %*% YY


  #RV index (Escoufier, 1973; Robert and Escoufier, 1976)
  RV <- tr(SXY %*% SYX)/sqrt(tr(SXX^2) * tr(SYY^2))
  #RLS index (Gower 1971; Lingoes and Schonemann (1974))
  RLS <- sqrt(tr(XX %*% t(XX) %*% YY %*% t(YY)))/ sqrt(tr(t(XX) %*% XX) %*% tr(t(YY) %*% YY))

  out <- data.frame(RV = RV, RLS = RLS)

  return(out)
}
