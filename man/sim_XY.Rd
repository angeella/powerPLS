% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_XY.R
\name{sim_XY}
\alias{sim_XY}
\title{simulate X data matrix under the alternative hypothesis}
\usage{
sim_XY(X, Y, A, nperm, A, scaling, post.transformation)
}
\arguments{
\item{out}{output from PLS2c}

\item{n}{number of observations}

\item{seed}{fix seed}

\item{post.transformation}{TRUE if you want to perform post transformation.}

\item{A}{number of components}
}
\value{
Returns a simulated X matrix under the alternative hypothesis.
}
\description{
simulate X data matrix under the alternative hypothesis.
}
\author{
Angela Andreella
}
