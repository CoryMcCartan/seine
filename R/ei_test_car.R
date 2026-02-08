#' Test the coarsening at random (CAR) assumption
#'
#' The CAR assumption, which is required for accurate estimates with [ei_est()],
#' implies that the conditional expectation function (CEF) of the aggregate
#' outcomes will take a certain partially linear form. This function tests that
#' implication by comparing a fully nonparametric estimate of the CEF to one
#' in the partially linear form implied by CAR, and comparing their goodness-of-fit.
#' It is **very important** that `preproc` be used in the `spec` to perform a
#' rich basis transformation of the covariates and predictors; a missing
#' `preproc` will lead to a warning.
#'
#' The test is carried out by fitting a regression on a fully basis-expanded
#' combination of covariates and predictors, and calculating the F statistic
#' compared to the "reduced" model that would be fit by [ei_ridge()]. To
#' account for penalization, the null distribution is estimated by permuting
#' the residuals from the reduced model.
#'
#' @param spec An `ei_spec` object created with [ei_spec()].
#'
#' @returns A 1-row tibble with columns describing the test results. The
#    `p.value` column contains the p-value for the test.
#'
#' @export
ei_test_car <- function(spec) {
}