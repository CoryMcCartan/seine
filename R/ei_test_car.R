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
#' by [ei_ridge()]. The test statistic is asymptotically chi-squared under the
#' null and may be anti-conservative in small samples, especially when the
#' dimensionality of the basis expansion is large.
#'
#' @param spec An `ei_spec` object created with [ei_spec()]. The object
#'   should use the `preproc` argument to specify a rich basis expansion of the
#'   covariates and predictors.
#' @inheritParams ei_ridge
#' @param iter The number of permutations to use when estimating the null
#'   distribution, including the original identity permutation.
#'   Ignored when `use_chisq = TRUE`.
#' @param undersmooth A value to divide the estimated ridge penalty by when
#'   partialling out the partially linear component of the model. A larger
#'   value will smooth the partially linear component less, which may improve
#'   Type I error control in finite samples at the cost of power.
#' @param use_chisq If `TRUE`, use the asymptotic chi-squared distribution for
#'   the Wald test statistic instead of conducting a permutation test. Only
#'   appropriate for larger sample sizes (Helwig 2022 recommends at least 200
#'   when a single predictor is used).
#' @param use_hc If `TRUE`, use a heteroskedasticity-consistent covariance estimate.
#'   More computationally intensive, but may make a difference in small samples
#'   or when there is substantial heteroskedasticity.
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
#' ei_test_car(spec, iter=20) # use a larger number in practice
#'
#' @export
ei_test_car <- function(spec, weights, iter = 1000, undersmooth = 1.5, use_chisq = nrow(spec) >= 2000, use_hc = FALSE) {
    validate_ei_spec(spec)
    if (!has_preproc(spec)) {
        cli_warn(c(
            "{.arg preproc} was not specified in your {.cls ei_spec} object",
            "i"="The {.fn ei_test_car} function relies on a rich basis expansion of the covariates and predictors.",
            "x"="Without a basis expansion in {.arg preproc}, the test will not be able to detect violations of the CAR assumption.",
            ">"="Consider basis expansions from the {.pkg bases} or {.pkg splines} package."
        ))
    }

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
    pen = fit0$penalty / undersmooth
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

    # observed test stat value
    fit = ridge_svd(udv, res, sqrt_w, pen, vcov = TRUE)

    # set up Wald stat calculation
    if (isTRUE(use_hc)) {
        xzr = (xzf - H0 %*% xzf)
        Sig_inv = solve((crossprod(xzr) + diag(ncol(xzr))*pen) / n)

        calc_wald = function(fit) {
            W = numeric(n_y)
            df = numeric(n_y)
            for (j in seq_len(n_y)) {
                Omega = crossprod((res - fit$fitted)[, j] * xzr) / n
                inv_vcov2 = pinv_sym(Sig_inv %*% Omega %*% Sig_inv)
                df[j] = attr(inv_vcov2, "rank")
                W[j] = n * crossprod(fit$coef[, j], inv_vcov2 %*% fit$coef[, j])
            }
            attr(W, "df") = df
            W
        }

        W0 = calc_wald(fit)
        df = attr(W0, "df")
    } else {
        inv_vcov = pinv_sym(fit$vcov_u)
        df = rep(attr(inv_vcov, "rank"), n_y)

        calc_wald = function(fit) {
            colSums(fit$coef * (inv_vcov %*% fit$coef)) / fit$sigma2
        }

        W0 = calc_wald(fit)
    }

    if (!isTRUE(use_chisq)) {
        if (!is.numeric(iter) && iter > 1) {
            cli_abort("{.arg iter} must be a positive integer.")
        }

        W = matrix(nrow = ncol(y), ncol = iter)
        W[, 1] = W0
        pb = cli::cli_progress_bar("Running permutations", total = iter)
        for (i in seq(2, iter, 1)) {
            res_p = res[sample.int(n), , drop=FALSE]
            fit_p = ridge_svd(udv, res_p, sqrt_w, pen, vcov = FALSE)
            W[, i] = calc_wald(fit_p)
            cli::cli_progress_update(id = pb)
        }
        cli::cli_progress_done(id = pb)
        p_val = rowSums(W >= W0) / iter
    } else {
        p_val = pchisq(W0, df = df, lower.tail = FALSE)
    }

    tibble::new_tibble(list(
        outcome = colnames(y),
        W = W0,
        df = df,
        p.value = p_val
    ))
}