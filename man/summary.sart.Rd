% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sart.R
\name{summary.sart}
\alias{summary.sart}
\alias{print.summary.sart}
\alias{print.sart}
\title{Summary for the Estimation of Tobit Models with Social Interactions}
\usage{
\method{summary}{sart}(object, Glist, data, ...)

\method{print}{summary.sart}(x, ...)

\method{print}{sart}(x, ...)
}
\arguments{
\item{object}{an object of class \code{sart}, output of the function \code{\link{sart}}.}

\item{Glist}{adjacency matrix or list sub-adjacency matrix. This is not necessary if the covariance method was computed in \link{cdnet}.}

\item{data}{dataframe containing the explanatory variables. This is not necessary if the covariance method was computed in \link{cdnet}.}

\item{...}{further arguments passed to or from other methods.}

\item{x}{an object of class \code{summary.sart}, output of the function \code{\link{summary.sart}}
or class \code{sart}, output of the function \code{\link{sart}}.}
}
\value{
A list of the same objects in \code{object}.
}
\description{
Summary and print methods for the class \code{sart} as returned by the function \link{sart}.
}
