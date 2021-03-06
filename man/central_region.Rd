% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/central_regions.R
\name{central_region}
\alias{central_region}
\title{Find the upper and lower bound for the \code{1-alpha} central region for a set of functions \code{fmat}}
\usage{
central_region(fmat, edepths, alpha = 0.05)
}
\arguments{
\item{fmat}{Matrix of functions. Each column is a function.}

\item{edepths}{Vector of extremal depths for the functions in fmat. first value should correspond
to the first function etc.}

\item{alpha}{Numeric. Default is 0.05.}
}
\value{
A list containing the upper and lower bound as separate vectors.
}
\description{
Find the upper and lower bound for the \code{1-alpha} central region for a set of functions \code{fmat}
}
