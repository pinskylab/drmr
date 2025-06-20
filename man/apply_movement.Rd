% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-gen.R
\name{apply_movement}
\alias{apply_movement}
\title{Simple movement to population dynamics}
\usage{
apply_movement(lambda, M, mov_age)
}
\arguments{
\item{lambda}{array of number of individuals per age, time, and patch}

\item{M}{movement matrix}

\item{mov_age}{ages at which movement starts (this can be generalized)}
}
\value{
an array of numbers by age, year and patch
}
\description{
Apply "simple" movement to population dynamics.
}
\author{
lcgodoy
}
