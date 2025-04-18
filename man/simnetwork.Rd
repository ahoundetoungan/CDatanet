% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxfunc.R
\name{simnetwork}
\alias{simnetwork}
\title{Simulating Network Data}
\usage{
simnetwork(dnetwork, normalise = FALSE)
}
\arguments{
\item{dnetwork}{A list of sub-network matrices, where the (i, j)-th position of the m-th matrix represents the probability that individual \code{i} is connected to individual \code{j} in the m-th network.}

\item{normalise}{A boolean indicating whether the returned matrices should be row-normalized (\code{TRUE}) or not (\code{FALSE}).}
}
\value{
A list of (row-normalized) adjacency matrices.
}
\description{
\code{simnetwork} generates adjacency matrices based on specified probabilities.
}
\examples{
# Generate a list of adjacency matrices
## Sub-network sizes
N         <- c(250, 370, 120)  
## Probability distributions
dnetwork  <- lapply(N, function(x) matrix(runif(x^2), x))
## Generate networks
G         <- simnetwork(dnetwork)
}
