# powerPLS
 
This package provides a tool to perform **power analysis** in Partial Least Squares (PLS) for classification when two classes are analyzed. 

## Installation

You can install the released version of `powerPLS` with:

``` r
devtools::install_github("angeella/powerPLS")
```

## Simulation

The main function is `computeSampleSize()` which estimated the power considering several values of sample size and number of latent components.

Here a toy example:

``` r
ind <- 0.1
ncomp <- 5
in.nvar <- 20  ##Number of variables
in.nclus <-2  ##Number of classes

## Within correlations
rho1 <- 0.2
rho2 <- 0.2
cormat <- list(group1 = diag(in.nvar) + rho1 - diag(in.nvar)*rho1,
               group2 = diag(in.nvar) + rho2 - diag(in.nvar)*rho2)

## Indicator validities

## Groups sizes for classes
in.clus.size <- c(15,15)

in.eta2 <- c(rep(ind, in.nvar/4), rep(0.001, in.nvar*3/4))
    
krt.lst <- lapply(seq(2), function(x) rnorm(in.nvar))
    
sim <- monte(seed = 123,
             nvar=in.nvar,
             nclus = in.nclus,
             clus.size = in.clus.size,
             eta2 = in.eta2,
             compactness = c(1, 1),
             kurt.list = krt.lst,
             cor.list = cormat)
```

``` r
X<- sim$data[, -1]
Y <- sim$data[,1]

computeSampleSize(X = X, Y = Y, A = ncomp,  n = 50,  alpha = 0.05, beta = 0.8,  = "auto-scaling",post.transformation = TRUE, Nsim = 100, nperm = 100, Y.prob = FALSE, test = "mcc",seed = 123)
```
