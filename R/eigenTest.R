#' @title randomization test
#' @description Performs randomization test based on eigenvalues
#' @usage eigenTest(X, Y, A, nperm, A, X.scaling, post.transformation)
#' @param X data matrix where columns represent the \eqn{p} classes and rows the \eqn{n} observations.
#' @param Y data matrix where columns represent the \eqn{k} variables and rows the \eqn{n} observations.
#' @param nperm number of permutations
#' @param A number of components
#' @param X.scaling
#' @param Y.prob
#' @param post.transformation TRUE if you want to apply post transformation.
#' @author Angela Andreella
#' @return Returns the corresponding pvalues
#' @importFrom compositions clrInv
#' @importFrom compositions clr
#' @export


eigenTest <- function(X, Y, nperm, A, X.scaling = "mean-centering", Y.prob = TRUE, eps = 0.01){

  if(X.scaling == "auto-scaling"){
    Mm <- apply(X, 2, mean)
    s <- apply(X, 2, sd)
    X <- (X - Mm)/s
  }
  if(X.scaling == "pareto-scaling"){
    Mm <- 0
    s <- apply(X, 2, sd)
    X <- (X - Mm)/s
  }
  if(X.scaling == "mean-centering"){
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
  E <- R <- R.p <- eigen.o <- w.p <- w <- r.p <- r <- Q <- list()
  pv <- c()
  eigen.p <- matrix(NA, ncol = nperm, nrow = A)

  E[[1]] <- as.matrix(X)
  R[[1]] <- as.matrix(Y)

  for(a in seq(A)){

    #Compute weight matrix
    out <- eigen(t(E[[a]]) %*% R[[a]] %*% t(R[[a]]) %*% E[[a]]) #well defined eigenvalue problem
    w[[a+1]] <- Re(out$vectors[,1])

   # w[[a+1]] <- (t(E[[a]]) %*% R[[a]])/(t(R[[a]] )%*% E[[a]] %*% t(E[[a]]) %*% R[[a]])[1]^(-1/2)
    r[[a+1]] <- E[[a]] %*% w[[a+1]]

    eigen.o[[a]] <- ((t(r[[a+1]]) %*% R[[a]])^2)[1]


    for(p in seq(nperm)){

      R.p[[p]] <- permuteIndex(R[[a]], by.row = TRUE, times = nrow(R[[a]]))
      out <- eigen(t(E[[a]]) %*% R.p[[p]] %*% t(R.p[[p]]) %*% E[[a]])
      w.p[[p]] <- Re(out$vectors[,1])
     # w.p[[p]] <- (t(E[[a]]) %*% R.p[[p]])/(t(R.p[[p]] )%*% E[[a]] %*% t(E[[a]]) %*%R.p[[p]])[1]^(-1/2)

      r.p[[p]] <- E[[a]] %*% w.p[[p]]

      eigen.p[a,p] <- ((t(r.p[[p]]) %*% R.p[[p]])^2)[1]

    }

    pv[a] <- (sum(eigen.p[a,] >= eigen.o[a]) + 1)/(nperm + 1)

    Q[[a+1]] <- diag(n) - r[[a+1]] %*% solve(t(r[[a+1]]) %*% r[[a+1]]) %*% t(r[[a+1]])
   # Q[[a+1]] <- diag(n) - r[[a+1]] %*% t(r[[a+1]])/ (t(r[[a+1]]) %*% r[[a+1]])[1]
    #Deflation step
    #X-deflation step
    E[[a+1]] <- Q[[a+1]] %*% E[[a]] #residual matrix X
    #Y-deflation step
    R[[a+1]] <- Q[[a+1]] %*% R[[a]] #residual matrix Y
  }

  if(is.unsorted(rev(unlist(eigen.o)))){#We will use only the predictive part to compute the eigenvalues

    warning("post-transformation is applied! Only the predictive part will be analyzed")
    #Rearrange weight matrix
    W <- NULL
    for (i in seq(length(w))) W <- cbind(W, w[[i]]) #number of obs times A
    out <- ptPLSc(X = X, Y = Y, W = W)
    Wtilde <- out$Wtilde
    M <- out$M
    T.score <- IDA(X = X, Y = Y, W = Wtilde)

    X.loading <- t(X) %*% T.score[,(M+1):A] %*% solve(t(T.score[,(M+1):A]) %*% T.score[,(M+1):A])

    X <- T.score[,(M+1):A] %*% t(X.loading)

    E[[1]] <- X
    R[[1]] <- Y

    for(a in seq(A)){

      #Compute weight matrix
      out <- eigen(t(E[[a]]) %*% R[[a]] %*% t(R[[a]]) %*% E[[a]]) #well defined eigenvalue problem
      w[[a+1]] <- Re(out$vectors[,1])

      # w[[a+1]] <- (t(E[[a]]) %*% R[[a]])/(t(R[[a]] )%*% E[[a]] %*% t(E[[a]]) %*% R[[a]])[1]^(-1/2)
      r[[a+1]] <- E[[a]] %*% w[[a+1]]

      eigen.o[[a]] <- ((t(r[[a+1]]) %*% R[[a]])^2)[1]


      for(p in seq(nperm)){

        R.p[[p]] <- permuteIndex(R[[a]], by.row = TRUE, times = nrow(R[[a]]))
        out <- eigen(t(E[[a]]) %*% R.p[[p]] %*% t(R.p[[p]]) %*% E[[a]])
        w.p[[p]] <- Re(out$vectors[,1])
        # w.p[[p]] <- (t(E[[a]]) %*% R.p[[p]])/(t(R.p[[p]] )%*% E[[a]] %*% t(E[[a]]) %*%R.p[[p]])[1]^(-1/2)

        r.p[[p]] <- E[[a]] %*% w.p[[p]]

        eigen.p[a,p] <- ((t(r.p[[p]]) %*% R.p[[p]])^2)[1]

      }

      pv[a] <- (sum(eigen.p[a,] >= eigen.o[a]) + 1)/(nperm + 1)

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
            pv_adj = pv_adj, test = unlist(eigen.o)))
}
