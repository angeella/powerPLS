% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeSampleSize.R
\name{computeSampleSize}
\alias{computeSampleSize}
\title{Sample size estimation}
\usage{
computeSampleSize(n, X, Y, A, alpha, beta,
nperm, Nsim, seed, test = "R2",...)
}
\arguments{
\item{n}{Vector of sample sizes to consider}

\item{X}{Data matrix where columns represent the \eqn{p} variables and
rows the \eqn{n} observations.}

\item{Y}{Data matrix where columns represent the two classes and
rows the \eqn{n} observations.}

\item{A}{Number of score components}

\item{alpha}{Type I error level. Default to 0.05}

\item{beta}{Type II error level. Default to 0.2.}

\item{nperm}{Number of permutations. Default to 100.}

\item{Nsim}{Number of simulations. Default to 100.}

\item{seed}{Seed value}

\item{test}{Type of test, one of \code{c("score", "mcc", "R2")}.
Default to "R2".}

\item{...}{Further parameters.}
}
\value{
Returns a data frame that contains the estimated power for each
sample size and number of components considered
}
\description{
Compute optimal sample size
}
\examples{
\dontrun{
datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
out <- computeSampleSize(X = datas$X, Y = datas$Y, A = 2, A = 3, n = 20, test = "R2")
}
}
\references{
For the general framework of power analysis for PLS-based methods see:

Andreella, A., Fino, L., Scarpa, B., & Stocchero, M. (2024). Towards a power analysis for PLS-based methods. arXiv preprint \url{https://arxiv.org/abs/2403.10289}.
}
\seealso{
\code{\link{computePower}}
}
\author{
Angela Andreella
}
