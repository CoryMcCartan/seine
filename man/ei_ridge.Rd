% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ei_ridge.R
\name{ei_ridge}
\alias{ei_ridge}
\alias{ei_ridge.formula}
\alias{ei_ridge.ei_spec}
\alias{ei_ridge.data.frame}
\alias{ei_ridge.matrix}
\alias{ei_ridge.default}
\title{Fit an ecological inference regression model}
\usage{
ei_ridge(x, ...)

\method{ei_ridge}{formula}(formula, data, weights, penalty = NULL, scale = TRUE, ...)

\method{ei_ridge}{ei_spec}(x, weights, penalty = NULL, scale = TRUE, ...)

\method{ei_ridge}{data.frame}(x, y, z, weights, penalty = NULL, scale = TRUE, ...)

\method{ei_ridge}{matrix}(x, y, z, weights, penalty = NULL, scale = TRUE, ...)

\method{ei_ridge}{default}(x, ...)
}
\arguments{
\item{x}{Depending on the context:
\itemize{
\item A \strong{data frame} of predictors.
\item A \strong{matrix} of predictors.
\item An \link{ei_spec} object containing the outcome, predictor, and covariates.
}

Predictors must be proportions that sum to 1 across rows.
You can use \code{\link[=ei_proportions]{ei_proportions()}} to assist in preparing predictor variables.
Covariates in an \link{ei_spec} object are shifted to have mean zero. If
\code{scale=TRUE} (the default), they are also scaled to have unit variance.}

\item{...}{Not currently used, but required for extensibility.}

\item{formula}{A formula such as \code{y ~ x0 + x1 | z} specifying the outcome \code{y}
regressed on the predictors of interest \code{x} and any covariates \code{z}.
The predictors should form a partition, that is, \code{x0 + x1 = 1} for each
observation. Users can be include more than two predictors as well, e.g.
\code{pct_white + pct_black + pct_hisp + pct_other}.
If there are just two predictors, it is acceptable to only include one in
the formula; the other will be formed as 1 minus the provided predictor.
Include additional covariates separated by a vertical bar \code{|}.
These covariates are strongly recommended for reliable ecological inference.
Covariates are shifted to have mean zero. If \code{scale=TRUE} (the default),
they are also scaled to have unit variance.}

\item{data}{When a \strong{formula} is used, \code{data} is a \strong{data frame} containing
both the predictors and the outcome.}

\item{weights}{<\code{\link[rlang:args_data_masking]{data-masking}}> A vector of unit
weights for estimation. These may be the same or different from the total
number of observations in each aggregate unit (see the \code{total} argument to
\code{\link[=ei_spec]{ei_spec()}}). See the discussion below under 'Weights' for choosing this
parameter. The default, uniform weights, makes a slightly
stronger-than-necessary assumption about the relationship between the
unit totals and the unknown data.}

\item{penalty}{The ridge penalty (a non-negative scalar). Set to \code{NULL} to
automatically determine the penalty which minimizes mean-square error,
via an efficient leave-one-out cross validation procedure.
The ridge regression solution is
\deqn{\hat\beta = (X^\top X + \lambda I)^{-1}X^\top y,}
where \eqn{\lambda} is the value of \code{penalty}.
Keep in mind when choosing \code{penalty} manually that covariates in \code{z} are
scaled to have mean zero and unit variance before fitting.}

\item{scale}{If \code{TRUE}, scale covariates \code{z} to have unit variance.}

\item{y}{When \code{x} is a \strong{data frame} or \strong{matrix}, \code{y} is the outcome
specified as:
\itemize{
\item A \strong{data frame} with numeric columns.
\item A \strong{matrix}
\item A numeric \strong{vector}.
}

When the outcome is a proportion, you can use \code{\link[=ei_proportions]{ei_proportions()}} to assist
in preparing it.}

\item{z}{When \code{x} is a \strong{data frame} or \strong{matrix}, \code{w} are any covariates,
specified as:
\itemize{
\item A \strong{data frame} with numeric columns.
\item A \strong{matrix}
}

These are shifted to have mean zero. If \code{scale=TRUE} (the default), they
are also scaled to have unit variance.}
}
\value{
An \code{ei_ridge} object, which supports various \link{ridge-methods}.
}
\description{
Fits a penalized regression model for ecological inference, allowing for
overall and unit-level estimates of conditional means using \code{\link[=ei_est]{ei_est()}}.
}
\details{
The regression is calculated using the singular value decomposition, which
allows for efficient recalculation under different \code{penalty} values as part
of leave-one-out cross-validation.
}
\section{Weights}{

The weakest identification result for ecological inference makes no
assumption about the number of observations per aggregate unit (the totals).
It requires, however, weighting the estimation steps according to the totals.
This may reduce efficiency when the totals are variable and a slightly
stronger condition holds.

Specifically, if the totals are conditionally mean-independent of the missing
data (the aggregation-unit level means of the outcome within each predictor
level), given covariates, then it is appropriate to use uniform weights in
estimation, or any fixed set of weights.

In general, estimation efficiency is improved when units with larger variance
in the outcome receive less weight. Various bulit-in options are provided by
the helper functions in \code{\link[=ei_wgt]{ei_wgt()}}.
}

\examples{
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(pop_urban, farm))
ei_ridge(spec)

ei_ridge(pres_dem_hum + pres_rep_nix + pres_ind_wal + pres_abs ~
      vap_white + vap_black + vap_other | pop_urban + farm, data = elec_1968)

}
