# powerPLS
 
This package provides a tool to perform **power analysis** in Partial Least Squares (PLS) for classification when two classes are analyzed. 

## Installation

You can install the released version of `powerPLS` with:

``` r
devtools::install_github("angeella/powerPLS")
```

## Quick overview

The main functions are 
- `computeSampleSize()` which estimated the power considering several values of sample size and number of score components.

- `computePower()` which estimated the power considering a fixed sample size and several number of score components.


``` r
datas <- simulatePilotData(nvar = 30, clus.size = c(5,5),m = 6,nvar_rel = 5,ncomp = 2)
out <- computePower(X = datas$X, Y = datas$Y, A = 3)
out <- computeSampleSize(X = datas$X, Y = datas$Y, A = 3, n = c(10,20,30))
```
