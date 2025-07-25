% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ei_local.R
\name{ei_est_local}
\alias{ei_est_local}
\alias{as.array.ei_est_local}
\title{Produce local ecological estimates}
\usage{
ei_est_local(
  regr,
  data,
  r_cov = NULL,
  bounds = NULL,
  conf_level = FALSE,
  unimodal = TRUE
)

\method{as.array}{ei_est_local}(x, ...)
}
\arguments{
\item{regr}{A fitted regression model, from \code{\link[=ei_ridge]{ei_ridge()}}, or another kind
of regression model wrapped with \code{\link[=ei_wrap_model]{ei_wrap_model()}}.}

\item{data}{The data frame, matrix, or \link{ei_spec} object that was used to fit
the regression.}

\item{r_cov}{A covariance matrix of the residuals to use in projecting the
local estimates onto the accounting constraint, or a list of matrices, one
for each outcome variable. Defaults to the identity matrix scaled by the
residual variance of \code{regr}, corresponding to orthogonal projection. Set
\code{r_cov=1} to use a degenerate covariance matrix corresponding to a (local)
neighborhood model. When there are multiple outcome variables and \code{r_cov} is
a matrix, it will be applied identically to each outcome.}

\item{bounds}{A vector \code{c(min, max)} of bounds for the outcome, to which the
local estimates will be truncated. In general, truncation will lead to
violations of the accounting identity. If \code{bounds = NULL}, they will be
inferred from the outcome variable: if it is contained within \eqn{[0, 1]},
for instance, then the bounds will be \code{c(0, 1)}. The default \code{bounds = FALSE} forces unbounded estimates.}

\item{conf_level}{A numeric specifying the level for confidence intervals.
If \code{FALSE} (the default), no confidence intervals are calculated.
For \code{regr} arguments from \code{\link[=ei_wrap_model]{ei_wrap_model()}}, confidence intervals will not
incorporate uncertainty in the prediction itself, just the residual. This
will trigger a warning periodically.}

\item{unimodal}{If \code{TRUE}, assume a unimodal residual distribution. Improves
width of confidence intervals by a factor of 4/9.}

\item{x}{An object of class \code{ei_est_local}}

\item{...}{Additional arguments (ignored)}
}
\value{
A data frame with estimates. The \code{.row} column in the output
corresponds to the observation index in the input. It has class
\code{ei_est_local}, supporting several methods.
}
\description{
Projects predictions from a fitted regression model onto the accounting
constraint using a provided residual covariance matrix. This ensures that
each set of local estimates satisfies the accounting identity. Local
estimates may be truncated to variable bounds.
}
\details{
Local estimates are produced independently for each outcome variable.
Truncation to bounds, if used, will in general lead to estimates that do
not satisfy the accounting identity.
}
\section{Methods (by generic)}{
\itemize{
\item \code{as.array(ei_est_local)}: Format estimates an array with dimensions
\verb{<rows>*<predictors>*<outcomes>}. Does not work if the object has been sorted.

}}
\examples{
\dontrun{
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(state, pop_urban, farm))

m = ei_ridge(spec)

ei_est_local(m, spec, conf_level = 0.95)
suppressWarnings(ei_est_local(m, spec, bounds=c(0.01, 0.2)))
}
}
