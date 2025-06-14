% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ei_ridge.R, R/ei_riesz.R
\name{ei_ridge_impl}
\alias{ei_ridge_impl}
\alias{ei_riesz_impl}
\title{Low-level implementations of \code{ei_ridge()} and \code{ei_riesz()}}
\usage{
ei_ridge_impl(x, y, z, weights, penalty = NULL)

ei_riesz_impl(x, z, total, weights, penalty)
}
\arguments{
\item{x}{A matrix of predictors}

\item{y}{A vector of outcomes}

\item{z}{A matrix of covariates}

\item{weights}{A vector of estimation weights}

\item{penalty}{The ridge penalty, which must be specified. Recommended value
is the same penalty used in \code{\link[=ei_ridge]{ei_ridge()}}, which is stored in the \code{penalty}
entry of the fitted model object.}

\item{total}{A vector of total observations per unit.}
}
\value{
A list with model components.
}
\description{
No checks are performed on the inputs.
Use of \code{\link[=ei_ridge]{ei_ridge()}} and \code{\link[=ei_riesz]{ei_riesz()}} is strongly recommended unless many
regressions must be fit, e.g., within a tight loop.
Only works for a single outcome, i.e., \code{y} must be a vector, not a matrix.
}
