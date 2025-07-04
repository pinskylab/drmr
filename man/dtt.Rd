% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialization.R
\name{dtt}
\alias{dtt}
\title{Density of a truncated Student's t distribution}
\usage{
dtt(x, mu = 0, sd = 1, df = 3, range = c(0, Inf), log = FALSE)
}
\arguments{
\item{x}{A \code{numeric vector} at which to evaluate the density.}

\item{mu}{A \code{numeric} scalar representing the mean of the underlying t
distribution. Defaults to 0.}

\item{sd}{A \code{numeric} scalar representing the standard deviation of the
underlying t distribution. Defaults to 1.}

\item{df}{A \code{numeric} scalar representing the degrees of freedom of the
underlying t distribution. Defaults to 3.}

\item{range}{A \code{numeric} vector of length 2 specifying the lower and
upper truncation bounds. Defaults to \code{c(0, Inf)}, indicating
truncation from below at 0.}

\item{log}{\code{Logical}; if \code{TRUE}, the natural logarithm of the
density is returned.  Defaults to \code{FALSE}.}
}
\value{
A numeric vector of the same length as \code{x}, containing the
(log) density values of the truncated t distribution.
}
\description{
Calculates the density of a truncated Student's t distribution.
}
\seealso{
\code{\link[stats:TDist]{stats::dt()}}, \code{\link[stats:TDist]{stats::pt()}}
}
\author{
lcgodoy
}
