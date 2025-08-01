% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ei_ridge.R
\name{ridge-methods}
\alias{ridge-methods}
\alias{predict.ei_ridge}
\alias{fitted.ei_ridge}
\alias{residuals.ei_ridge}
\alias{vcov.ei_ridge}
\alias{summary.ei_ridge}
\alias{weights.ei_ridge}
\title{Methods for \code{\link{ei_ridge}} models}
\usage{
\method{predict}{ei_ridge}(object, new_data, type = "numeric", ...)

\method{fitted}{ei_ridge}(object, ...)

\method{residuals}{ei_ridge}(object, ...)

\method{vcov}{ei_ridge}(object, ...)

\method{summary}{ei_ridge}(object, ...)

\method{weights}{ei_ridge}(object, normalize = TRUE, ...)
}
\arguments{
\item{object}{A fitted \link{ei_ridge} model}

\item{new_data}{A data frame, matrix, or \link{ei_spec} of new predictors.}

\item{type}{The type of predictions to generate; only \code{"numeric"} is supported.}

\item{...}{Additional arguments (ignored)}

\item{normalize}{If \code{TRUE}, normalize the weights to have mean 1.}
}
\description{
Models fitted with \code{\link[=ei_ridge]{ei_ridge()}} support various generic methods.
}
\section{Functions}{
\itemize{
\item \code{predict(ei_ridge)}: Predict from an \code{ei_ridge} model.

\item \code{fitted(ei_ridge)}: Extract fitted values.

\item \code{residuals(ei_ridge)}: Extract residuals.

\item \code{vcov(ei_ridge)}: Extract covariance of coefficient estimates.

\item \code{summary(ei_ridge)}: Summarize the model's fitted values and \eqn{R^2}.

\item \code{weights(ei_ridge)}: Extract estimation weights from a fitted model.

}}
