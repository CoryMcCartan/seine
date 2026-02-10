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
#' The test is a Kennedy-Cade (1996) style permutation test on a Wald statistic
#' for the coefficients not included in the "reduced" model that would be fit
#' by [ei_ridge()].
#' The test is carried out by fitting a regression on a fully basis-expanded
#' combination of covariates and predictors, and calculating a Wald statistic
#' for the
#'
#' @param spec An `ei_spec` object created with [ei_spec()]. The object
#'   should use the `preproc` argument to specify a rich basis expansion of the
#'   covariates and predictors.
#' @inheritParams ei_ridge
#' @param iter The number of permutations to use when estimating the null
#'   distribution. Ignored when `use_chisq = TRUE`.
#' @param use_chisq If `TRUE`, use the asymptotic chi-squared distribution for
#'   the Wald test statistic instead of conducting a permutation test. Only
#'   appropriate for larger sample sizes (Helwig 2022 recommends at least 200
#'   when a single predictor is used).
#'
#' @returns A tibble with one row per outcome variable and columns describing
#'    the test results. The `p.value` column contains the p-values for the test.
#'    P-values are not adjusted by default; pass them to [stats::p.adjust()] if
#'    desired.
#'
#' @references
#' Helwig, N. E. (2022). Robust Permutation Tests for Penalized Splines. _Stats_,
#' 5(3), 916-933.
#'
#' Kennedy, P. E., & Cade, B. S. (1996). Randomization tests for multiple regression.
#' _Communications in Statistics-Simulation and Computation_, 25(4), 923-936.
#'
#' McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
#' estimation of conditional means from aggregate data.
#' Working paper [arXiv:2509.20194](https://arxiv.org/abs/2509.20194).
#'
#' @examples
#' data(elec_1968)
#'
#' # basis expansion: poly() with degree=2 not recommended in practice
#' preproc = if (requireNamespace("bases", quietly = TRUE)) {
#'     ~ bases::b_bart(.x, trees = 100)
#' } else {
#'     ~ poly(as.matrix(.x), degree=2, simple=TRUE)
#' }
#'
#' spec = ei_spec(
#'     data = elec_1968,
#'     predictors = vap_white:vap_other,
#'     outcome = pres_dem_hum:pres_abs,
#'     total = pres_total,
#'     covariates = c(pop_city:pop_rural, farm:educ_coll, starts_with("inc_")),
#'     preproc = preproc
#' )
#'
#' ei_test_car(spec, iter=19) # use a larger number in practice
#'
#' @export
ei_test_car <- function(spec, weights, iter = 1000, use_chisq = FALSE) {
    validate_ei_spec(spec)
    n = nrow(spec)
    x_col = attr(spec, "ei_x")
    z_col = attr(spec, "ei_z")
    x = spec[, x_col, drop = FALSE]
    z = spec[, z_col, drop = FALSE]
    z_proc = attr(spec, "ei_z_proc")

    int_scale = 1e5
    xz0 = row_kronecker(as.matrix(x), z_proc, int_scale)
    xzf = run_preproc(spec, z_col = c(x_col, z_col))

    if (missing(weights)) {
        weights = rep(1, n)
    } else {
        weights = eval_tidy(enquo(weights), spec)
    }
    sqrt_w = sqrt(weights / mean(weights))

    y = as.matrix(spec[, attr(spec, "ei_y"), drop = FALSE])
    n_y = ncol(y)

    # first, residualize out xz0
    udv0 = svd(xz0 * sqrt_w)
    fit0 = ridge_auto(udv0, y, sqrt_w, vcov = FALSE)
    pen = fit0$penalty
    d_pen_h = udv0$d^2 / (udv0$d^2 + pen)
    H0 = tcrossprod(scale_cols(udv0$u, d_pen_h), udv0$u)
    res = y - H0 %*% y
    udv = svd((xzf - H0 %*% xzf) * sqrt_w)

    # pseudo-inverse
    pinv_sym = function(M) {
        eig = eigen(M, symmetric=TRUE)
        rk = seq_len(sum(eig$values > 1e-10))
        eig$values[rk] = 1 / eig$values[rk]
        out = tcrossprod(scale_cols(eig$vectors, eig$values), eig$vectors)
        attr(out, "rank") = length(rk)
        out
    }

    fit = ridge_svd(udv, res, sqrt_w, pen, vcov = TRUE)
    inv_vcov = pinv_sym(fit$vcov_u)

    calc_wald = function(fit) {
        colSums((inv_vcov %*% fit$coef) * fit$coef) / fit$sigma2
    }
    W0 = calc_wald(fit)

    if (!isTRUE(use_chisq)) {
        if (!is.numeric(iter) && iter > 0) {
            cli_abort("{.arg iter} must be a positive integer.")
        }

        W = matrix(nrow = ncol(y), ncol = iter)
        pb = cli::cli_progress_bar("Running permutations", total = iter)
        for (i in seq_len(iter)) {
            res_p = res[sample.int(n), , drop=FALSE]
            fit_p = ridge_svd(udv, res_p, sqrt_w, pen, vcov = FALSE)
            W[, i] = calc_wald(fit_p)
            cli::cli_progress_update(id = pb)
        }
        cli::cli_progress_done(id = pb)
        p_val = (rowSums(W >= W0) + 1) / (iter + 1)
    } else {
        p_val = pchisq(W0, df = attr(inv_vcov, "rank"), lower.tail = FALSE)
    }

    tibble::new_tibble(list(
        outcome = colnames(y),
        W = W0,
        df = attr(inv_vcov, "rank"),
        p.value = p_val
    ))
}