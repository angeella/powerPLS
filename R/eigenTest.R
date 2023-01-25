#' @title randomization test
#' @description Performs randomization test based on eigenvalues
#' @usage eigenTest(X, Y, nperm, A, scaling = "mean-centering", Y.prob = TRUE, eps = 0.01)
#' @param X data matrix where columns represent the \eqn{p} classes and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and rows the \eqn{n} observations.
#' @param nperm number of permutations
#' @param A number of components
#' @param scaling type of scaling
#' @param Y.prob TRUE if Y describes the probability to being in the class.
#' @param eps see details
#' @author Angela Andreella
#' @return Returns the corresponding pvalues
#' @importFrom compositions clrInv
#' @importFrom compositions clr
#' @importFrom stats sd
#' @export


eigenTest <- function(X, Y, nperm, A, scaling = "mean-centering", Y.prob = TRUE, eps = 0.01){

  if(scaling == "auto-scaling"){
    Mm <- apply(X, 2, mean)
    s <- apply(X, 2, sd)
    X <- (X - Mm)/s
  }
  if(scaling == "pareto-scaling"){
    Mm <- 0
    s <- apply(X, 2, sd)
    X <- (X - Mm)/s
  }
  if(scaling == "mean-centering"){
    Mm <- apply(X, 2, mean)
    s <- 1
    X <- (X - Mm)/s
  }

  if(!Y.prob){
    Y[which(Y==0)]<-eps
    Y[which(Y==1)]<-1-(ncol(Y)-1)*eps
    P <- matrix(clr(Y), ncol = ncol(Y))
    Mm <- apply(P, 2, mean)
    s <- apply(P, 2, sd)
    P <- (P - Mm)/s
    Y <- as.matrix(P)
  }

  n <- nrow(Y)
  E <- R <- R_p <- eigen_o <- w_p <- w <- r_p <- r <- Q <- list()
  pv <- c()
  eigen_p <- matrix(NA, ncol = nperm, nrow = A)

  E[[1]] <- as.matrix(X)
  R[[1]] <- as.matrix(Y)

  for(a in seq(A)){

    #Compute weight matrix
    out <- eigen(t(E[[a]]) %*% R[[a]] %*% t(R[[a]]) %*% E[[a]]) #well defined eigenvalue problem
    w[[a+1]] <- Re(out$vectors[,1])

   # w[[a+1]] <- (t(E[[a]]) %*% R[[a]])/(t(R[[a]] )%*% E[[a]] %*% t(E[[a]]) %*% R[[a]])[1]^(-1/2)
    r[[a+1]] <- E[[a]] %*% w[[a+1]]

    eigen_o[[a]] <- ((t(r[[a+1]]) %*% R[[a]])^2)[1]


    for(p in seq(nperm)){

      R_p[[p]] <- permuteIndex(R[[a]], by.row = TRUE, times = nrow(R[[a]]))
      out <- eigen(t(E[[a]]) %*% R_p[[p]] %*% t(R_p[[p]]) %*% E[[a]])
      w_p[[p]] <- Re(out$vectors[,1])
     # w_p[[p]] <- (t(E[[a]]) %*% R_p[[p]])/(t(R_p[[p]] )%*% E[[a]] %*% t(E[[a]]) %*%R_p[[p]])[1]^(-1/2)

      r_p[[p]] <- E[[a]] %*% w_p[[p]]

      eigen_p[a,p] <- ((t(r_p[[p]]) %*% R_p[[p]])^2)[1]

    }

    pv[a] <- (sum(eigen_p[a,] >= eigen_o[a]) + 1)/(nperm + 1)

    Q[[a+1]] <- diag(n) - r[[a+1]] %*% solve(t(r[[a+1]]) %*% r[[a+1]]) %*% t(r[[a+1]])
   # Q[[a+1]] <- diag(n) - r[[a+1]] %*% t(r[[a+1]])/ (t(r[[a+1]]) %*% r[[a+1]])[1]
    #Deflation step
    #X-deflation step
    E[[a+1]] <- Q[[a+1]] %*% E[[a]] #residual matrix X
    #Y-deflation step
    R[[a+1]] <- Q[[a+1]] %*% R[[a]] #residual matrix Y
  }

  if(is.unsorted(rev(unlist(eigen_o)))){#We will use only the predictive part to compute the eigenvalues

    warning("post-transformation is applied! Only the predictive part will be analyzed")
    #Rearrange weight matrix
    W <- NULL
    for (i in seq(length(w))) W <- cbind(W, w[[i]]) #number of obs times A
    out <- ptPLSc(X = X, Y = Y, W = W)
    Wtilde <- out$Wtilde
    M <- out$M
    T_score <- IDA(X = X, Y = Y, W = Wtilde)

    if(A==M){
      X_loading <- t(X) %*% T_score %*% solve(t(T_score) %*% T_score)
      X <- T_score %*% t(X_loading)
    }else{
      X_loading <- t(X) %*% T_score[,(M+1):A] %*% solve(t(T_score[,(M+1):A]) %*% T_score[,(M+1):A])
      X <- T_score[,(M+1):A] %*% t(X_loading)
    }


    E[[1]] <- X
    R[[1]] <- Y

    for(a in seq(A)){

      #Compute weight matrix
      out <- eigen(t(E[[a]]) %*% R[[a]] %*% t(R[[a]]) %*% E[[a]]) #well defined eigenvalue problem
      w[[a+1]] <- Re(out$vectors[,1])

      # w[[a+1]] <- (t(E[[a]]) %*% R[[a]])/(t(R[[a]] )%*% E[[a]] %*% t(E[[a]]) %*% R[[a]])[1]^(-1/2)
      r[[a+1]] <- E[[a]] %*% w[[a+1]]

      eigen_o[[a]] <- ((t(r[[a+1]]) %*% R[[a]])^2)[1]


      for(p in seq(nperm)){

        R_p[[p]] <- permuteIndex(R[[a]], by.row = TRUE, times = nrow(R[[a]]))
        out <- eigen(t(E[[a]]) %*% R_p[[p]] %*% t(R_p[[p]]) %*% E[[a]])
        w_p[[p]] <- Re(out$vectors[,1])
        # w_p[[p]] <- (t(E[[a]]) %*% R_p[[p]])/(t(R_p[[p]] )%*% E[[a]] %*% t(E[[a]]) %*%R_p[[p]])[1]^(-1/2)

        r_p[[p]] <- E[[a]] %*% w_p[[p]]

        eigen_p[a,p] <- ((t(r_p[[p]]) %*% R_p[[p]])^2)[1]

      }

      pv[a] <- (sum(eigen_p[a,] >= eigen_o[a]) + 1)/(nperm + 1)

      Q[[a+1]] <- diag(n) - r[[a+1]] %*% solve(t(r[[a+1]]) %*% r[[a+1]]) %*% t(r[[a+1]])
      # Q[[a+1]] <- diag(n) - r[[a+1]] %*% t(r[[a+1]])/ (t(r[[a+1]]) %*% r[[a+1]])[1]
      #Deflation step
      #X-deflation step
      E[[a+1]] <- Q[[a+1]] %*% E[[a]] #residual matrix X
      #Y-deflation step
      R[[a+1]] <- Q[[a+1]] %*% R[[a]] #residual matrix Y
    }

  }
  #Adjust p-values
  pv_adj <- sapply(c(1:A), function(x) max(pv[1:x]))

return(list(pv = pv,
            pv_adj = pv_adj, test = unlist(eigen_o)))
}
