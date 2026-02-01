#' Produce local ecological estimates
#'
#' Projects predictions from a fitted regression model onto the accounting
#' constraint using a provided residual covariance matrix. This ensures that
#' each set of local estimates satisfies the accounting identity. Local
#' estimates may be truncated to variable bounds.
#'
#' Local estimates are produced jointly across outcome variables. When bounds
#' are applied, unless `sum_one = TRUE`, the estimates for each observation may
#' not satisfy logical constraints, including the accounting identity.
#'
#' Projections are done obliquely in accordance with `r_cov` via quadratic
#' programming. Occasionally, the quadratic program may be infeasible due to
#' the specific data, features of `r_cov`, or numerical errors. Various
#' relaxations of the accounting identity and `r_cov` are attempted in these cases;
#' indices where relaxations of `r_cov` were used are stored in the `proj_relax`
#' attribute of the output, and indices of infeasible projections are stored in
#' the `proj_misses` attribute.
#'
#' @param regr A fitted regression model, from [ei_ridge()], or another kind
#'    of regression model wrapped with [ei_wrap_model()].
#' @param data The data frame, matrix, or [ei_spec] object that was used to fit
#'   the regression.
#' @param total <[`tidy-select`][tidyselect::select_helpers]> A variable
#'   containing the total number of observations in each aggregate unit. For
#'   example, the column containing the total number of voters. Required if
#'   `data` is not an [ei_spec()] object.
#' @param r_cov A covariance matrix of the residuals to use in projecting the
#'   local estimates onto the accounting constraint, such as one estimated with
#'   [ei_resid_cov()]. Defaults to the identity matrix scaled by the residual
#'   variance of `regr`, corresponding to orthogonal projection.
#'   Set `r_cov=1` to use a degenerate covariance matrix corresponding to a
#'   (local) neighborhood model. When there are multiple outcome variables and
#'   r_cov` is a matrix with entries for each predictor, it will be applied
#'   identically to each outcome. Alternatively, a matrix with entries for each
#'   predictor-outcome combination may be provided, with entries in the order
#'   (Y1|X1, Y1|X2, ..., Y2|X1, Y2|X2, ...).
#' @param bounds A vector `c(min, max)` of bounds for the outcome, to which the
#'   local estimates will be truncated. In general, truncation will lead to
#'   violations of the accounting identity. If `bounds = NULL`, they will be
#'   inferred from the outcome variable: if it is contained within \eqn{[0, 1]},
#'   for instance, then the bounds will be `c(0, 1)`. Setting `bounds = FALSE`
#'   forces unbounded estimates. The default uses the `bounds` attribute of
#'   `regr`, if available, or infers from the outcome variable otherwise.
#'   Note that bounds are currently not applied if `contrast` is provided.
#' @inheritParams ei_ridge
#' @inheritParams ei_est
#' @param conf_level A numeric specifying the level for confidence intervals.
#'   If `FALSE` (the default), no confidence intervals are calculated.
#'   For `regr` arguments from [ei_wrap_model()], confidence intervals will not
#'   incorporate uncertainty in the prediction itself, just the residual. This
#'   will trigger a warning periodically.
#' @param regr_var If `TRUE`, incorporate uncertainty from the regression model
#'   when calculating confidence intervals. Only applies when `regr` is fitted
#'   with [ei_ridge()], and requires that function be called with `vcov = TRUE`.
#' @param unimodal If `TRUE`, assume a unimodal residual distribution. Improves
#'   width of confidence intervals by a factor of 4/9.
#'
#' @returns A data frame with estimates. The `.row` column in the output
#'   corresponds to the observation index in the input. The `wt` column contains
#'   the product of the predictor variable and total for each observation.
#'   Taking a weighted average of the estimate against this column will produce
#'   a global estimate. It has class `ei_est_local`, supporting several methods.
#'
#' @inherit ei_est references
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
#'                total = pres_total, covariates = c(state, pop_urban, farm))
#'
#' m = ei_ridge(spec)
#'
#' ei_est_local(m, spec, bounds = c(0, 1), sum_one = TRUE, conf_level = 0.95)
#'
#' r_cov = ei_resid_cov(m, spec)
#' e_orth = ei_est_local(m, spec, bounds = c(0, 1), sum_one = TRUE, conf_level = 0.95)
#' e_nbhd = ei_est_local(m, spec, r_cov = 1, bounds = c(0, 1), sum_one = TRUE, conf_level = 0.95)
#' e_rcov = ei_est_local(m, spec, r_cov = r_cov, bounds = c(0, 1), sum_one = TRUE, conf_level = 0.95)
#' # average interval width
#' c(
#'     e_orth = mean(e_orth$conf.high - e_orth$conf.low),
#'     e_nbhd = mean(e_nbhd$conf.high - e_nbhd$conf.low),
#'     e_rcov = mean(e_rcov$conf.high - e_rcov$conf.low)
#' )
#' @export
ei_est_local = function(
    regr,
    data,
    total,
    r_cov = NULL,
    contrast = NULL,
    bounds = regr$blueprint$bounds,
    sum_one = NULL,
    conf_level = FALSE,
    regr_var = TRUE,
    unimodal = TRUE
) {
    y = est_check_outcome(regr, data, NULL)
    n = nrow(y)
    n_y = ncol(y)

    cli_warn(
        "Local confidence intervals do not yet incorporate prediction uncertainty.",
        .frequency = "regularly",
        .frequency_id = "ei_est_local_temp"
    )
    if (inherits(regr, "ei_wrapped") && isTRUE(regr_var) && !isFALSE(conf_level)) {
        cli_warn(
            "Local confidence intervals with wrapped model objects
                  do not incorporate prediction uncertainty.",
            .frequency = "regularly",
            .frequency_id = "ei_est_local"
        )
        regr_var = FALSE
    }

    rl = est_check_regr(regr, data, n, NULL, n_y, vcov = isTRUE(regr_var))
    n_x = length(rl$preds)

    bounds = ei_bounds(bounds, y, clamp = 1e-8)
    if (is.null(sum_one) && all(bounds == c(0, 1))) {
        sum_one = isTRUE(all.equal(rowSums(y), rep(1, nrow(y))))
    }

    if (missing(total)) {
        if (inherits(data, "ei_spec")) {
            total = attr(data, "ei_n")
        } else {
            cli_abort("{.arg total} is required when {.arg data} is not an {.cls ei_spec} object.")
        }
    } else {
        total = eval_tidy(enquo(total), data)
    }

    r_cov = check_proc_r_cov(r_cov, regr$sigma2, n_x)

    contr = check_contrast(contrast, colnames(rl$x), colnames(y))
    x_nm = contr$x_nm
    y_nm = contr$y_nm

    eta = matrix(nrow = n, ncol = n_x * n_y)
    for (k in seq_len(n_x)) {
        idx = k + seq(0, by = n_x, length.out = n_y)
        eta[, idx] = rl$preds[[k]]
    }
    eta_proj = local_proj(rl$x, eta, y - rl$yhat, r_cov, bounds, sum_one)
    if (!is.null(contrast)) {
        eta_proj = eta_proj %*% contr$m
    }

    sds = if (!isFALSE(conf_level)) {
        local_sds(rl$x, r_cov, rl$vcov, contr$m, !is.null(contrast))
    } else {
        NULL
    }

    ests = list(
        .row = rep(seq_len(n), length(x_nm)),
        predictor = rep(x_nm, each = n),
        outcome = rep(y_nm, each = n),
        wt = if (is.null(contrast)) rep(c(rl$x * total), n_y) else NULL,
        estimate = c(eta_proj),
        std.error = c(sds)
    )
    if (!is.null(contrast)) {
        ests[["wt"]] = NULL
    }
    ests = tibble::new_tibble(
        ests,
        proj_misses = attr(eta_proj, "misses"),
        proj_relax = attr(eta_proj, "relax"),
        class = "ei_est_local"
    )

    if (!isFALSE(conf_level)) {
        fac = if (isTRUE(unimodal)) 4 / 9 else 1
        chebyshev = sqrt(fac / (1 - conf_level))
        ests$conf.low = ests$estimate - chebyshev * ests$std.error
        ests$conf.high = ests$estimate + chebyshev * ests$std.error

        if (is.null(contrast)) {
            ests$conf.low[ests$conf.low < bounds[1]] = bounds[1]
            ests$conf.high[ests$conf.high < bounds[1]] = bounds[1]
            ests$conf.low[ests$conf.low > bounds[2]] = bounds[2]
            ests$conf.high[ests$conf.high > bounds[2]] = bounds[2]
        }
    }
    ests$std.error = NULL

    ests
}

#' @describeIn ei_est_local Format estimates an array with dimensions
#'   `<rows>*<predictors>*<outcomes>`. Does not work if the object has been sorted.
#' @param x An object of class `ei_est_local`
#' @param ... Additional arguments (ignored)
#' @export
as.array.ei_est_local = function(x, ...) {
    nm_x = unique(x$predictor)
    nm_y = unique(x$outcome)
    n_x = length(nm_x)
    n_y = length(nm_y)
    n = nrow(x) / n_x / n_y
    array(x$estimate, dim=c(n, n_x, n_y), dimnames=list(NULL, nm_x, nm_y))
}

#' Estimate the residual covariance of the local estimands
#'
#' Under a slightly stronger coarsening at random assumption (applying to
#' second moments), _and_ an assumption of homoskedasticity in the covariates,
#' this function estimates the covariance matrix of the local estimands
#' \eqn{\beta_{gj}=\mathbb{E}[Y | X_j=1, Z=z_g, G=g]} around their local mean.
#' See the reference for more detail.
#'
#' Homoskedasticity in the covariates implies that the variance of the residuals
#' depends linearly on the entries of \eqn{\bar{X}\bar{X}^\top}. This function
#' fits an auto-tuned ridge regression of the empirical second moments of the
#' residuals on these predictors, and uses the polarization identity discussed
#' in the references to estimate the covariance for each local estimand. When
#' the estiamated covariance is not positive semidefinite, it is projected onto
#' the cone of positive semidefinite matrices.
#'
#' @param regr A fitted regression model, from [ei_ridge()], or another kind
#'    of regression model wrapped with [ei_wrap_model()].
#' @param data The data frame, matrix, or [ei_spec] object that was used to fit
#'   the regression.
#'
#' @returns A covariance matrix. The variables are ordered by predictor within
#'   outcome, e.g. (Y1|X1, Y1|X2, ..., Y2|X1, Y2|X2, ...).
#'
#' @inherit ei_est references
#'
#' @export
ei_resid_cov <- function(regr, data) {
    if (!inherits(regr, "ei_ridge")) {
        cli_abort("{.fun ei_resid_cov} only supports regressions fit with {.fn ei_ridge}.")
    }
    y = regr$y
    n = nrow(y)
    n_y = ncol(y)

    xcols = regr$blueprint$ei_x
    idx_x = match(xcols, colnames(data))
    x = as.matrix(pull_x(data, idx_x))
    n_x = length(xcols)

    idx_tri = c(lower.tri(diag(n_y), diag = TRUE))
    yr = row_kronecker(resid(regr), resid(regr), 0)[, idx_tri]
    udv = svd(cbind(row_kronecker(x, x, 0)))

    fit = ridge_auto(udv, yr, rep(1, n), vcov = FALSE)
    pred_sigma = function(j, k) {
        x_plug = numeric(n_x)
        x_plug[c(j, k)] = 1
        sigma = matrix(0, n_y, n_y)
        pred = 0.5 * (x_plug %x% x_plug) %*% fit$coef
        sigma[lower.tri(sigma, diag = TRUE)] = pred
        sigma[upper.tri(sigma)] = t(sigma)[upper.tri(sigma)]
        sigma
    }

    r_cov0 = matrix(0, n_x * n_y, n_x * n_y)
    for (j in seq_len(n_x)) {
        for (k in seq_len(j)) {
            off = (seq_len(n_y) - 1) * n_x
            pred_part = 2 * (pred_sigma(j, k) - 0.25 * pred_sigma(j, j) - 0.25 * pred_sigma(k, k))
            r_cov0[j + off, k + off] = pred_part
            r_cov0[k + off, j + off] = pred_part
        }
    }
    eig = eigen(r_cov0)
    r_cov = eig$vectors %*% diag(pmax(eig$values, 1e-8)) %*% t(eig$vectors)

    cov_nm = c(outer(colnames(x), colnames(y), paste, sep=":"))
    colnames(r_cov) <- rownames(r_cov) <- cov_nm

    r_cov
}

# Process r_cov
check_proc_r_cov <- function(r_cov, sigma2, n_x) {
    n_y = length(sigma2)
    if (is.matrix(r_cov)) {
        if (all(dim(r_cov) == n_x * n_y)) {
            # ok
        } else if (all(dim(r_cov) == n_x)) {
            r_cov = diag(n_y) %x% r_cov
        } else {
            cli_abort(c(
                "Invalid {.arg r_cov} matrix dimensions.",
                ">" = "Expected either {n_x}x{n_x} or {n_x * n_y}x{n_x * n_y}."
            ), call = parent.frame())
        }
    } else if (is.null(r_cov)) {
        r_cov = diag(sigma2) %x% diag(n_x)
    } else if (length(r_cov) == 1 && r_cov == 1) {
        r_cov = diag(sigma2) %x% (1 + diag(n_x) * 1e-8)
    } else {
        cli_abort("Invalid {.arg r_cov} format. Consult documentation.", call = parent.frame())
    }
    chol(r_cov)
}

# Solve QP to project estimates onto tomography plane and into bounds
# Not the fastest possible implementation (pure C++ would be better), but fast enough
# r_cov here is Cholesky factor
local_proj = function(x, eta, eps, r_cov, bounds, sum_one) {
    n = nrow(eta)
    n_x = ncol(x)
    n_y = ncol(eps)
    sum_one = isTRUE(sum_one)
    eta_diff = matrix(nrow = n, ncol = n_x * n_y)

    # avoid overflow
    r_cov = r_cov / sqrt(norm(crossprod(r_cov), "2"))

    # parameters are the displacement in each estimate
    # (x1y1, x2y1, x3y1, x1y2, x2y2, x3y2, ...)
    # minimize overall displacement st x-weighted displacement = residual
    # and (optionally) bounds and sum-to-one constraints are satisfied
    zeros = rep(0, n_x * n_y)
    Amat = matrix(0, nrow = n_x * n_y, ncol = n_y * 2) # i-specific, filled later
    rownames(Amat) = colnames(eta) # helps debug
    b0 = cbind(eps, -eps)
    if (sum_one) {
        if (n_y == 1 || all(bounds == c(-Inf, Inf))) {
            cli_abort(
                "Using{.arg sum_one} requires multiple bounded outcomes.",
                call = parent.frame()
            )
        }
        rs_mat =  rep(1, n_y) %x% diag(n_x)
        Amat = cbind(rs_mat, Amat)
        b0 = cbind(1 - eta %*% rs_mat, b0)
    }
    if (!is.infinite(bounds[1])) {
        Amat = cbind(Amat, diag(n_x * n_y))
        b0 = cbind(b0, bounds[1] - eta)
    }
    if (!is.infinite(bounds[2])) {
        Amat = cbind(Amat, -diag(n_x * n_y))
        b0 = cbind(b0, -bounds[2] + eta)
    }

    idx_eps = sum_one * n_x + seq_len(2*n_y)
    patt_eps = cbind(diag(n_y), -diag(n_y))

    constr_pt = function(Dmat, bvec, tol) {
        bvec[idx_eps] = bvec[idx_eps] - tol
        quadprog::solve.QP(
            Dmat = Dmat, # distance metric
            dvec = zeros,
            Amat = Amat,
            bvec = bvec,
            meq = sum_one * n_x,
            factorized = TRUE
        )$solution
    }

    misses = integer(0)
    relaxations = 0
    r_cov_relax = diag(diag(r_cov) + 1e-3)
    for (i in seq_len(n)) {
        Amat[, idx_eps] = patt_eps %x% x[i, ]
        tol = 1e-12
        relax_D = FALSE
        repeat {
            Dmat = if (!relax_D) r_cov else r_cov_relax
            ans = tryCatch(constr_pt(Dmat, b0[i, ], tol), error = \(e) NULL)
            if (!is.null(ans)) break
            if (tol > 0.0005) {
                if (!relax_D) {
                    relaxations = relaxations + 1
                    relax_D = TRUE
                    tol = 1e-12
                    next
                }
                misses <- c(misses, i)
                ans = rep(eps[i, ], n_x)
                break
            }
            tol = tol * 1000
        }
        eta_diff[i, ] = ans
    }

    out = eta + eta_diff
    if (!is.infinite(bounds[1])) {
        out[out < bounds[1]] = bounds[1]
    }
    if (!is.infinite(bounds[2])) {
        out[out > bounds[2]] = bounds[2]
    }
    attr(out, "relax") = relaxations
    attr(out, "misses") = misses
    out
}

# r_cov here is Cholesky factor
local_sds = function(x, r_cov, regr_cov, contr_m, is_contr = FALSE) {
    n = nrow(x)
    n_x = ncol(x)
    n_y = nrow(r_cov) / n_x
    sds = matrix(nrow = n, ncol = ncol(contr_m))
    for (i in seq_len(n)) {
        H = diag(n_y) %x% local_basis(x[i, ])
        R_i = if (!is.null(regr_cov)) {
            chol(crossprod(r_cov) + (diag(n_y) %x% matrix(regr_cov[i, ], n_x, n_x)))
        } else {
            r_cov
        }
        Pi_Sigma = oblique_proj(H, R_i)

        if (is_contr) {
            sds[i, ] = sqrt(colSums((Pi_Sigma %*% contr_m) * contr_m))
        } else {
            sds[i, ] = sqrt(diag(Pi_Sigma))
        }
    }
    sds
}

local_basis = function(x) {
    qr.Q(qr(cbind(x, diag(length(x)))))[, -1, drop=FALSE]
}

# project onto plane spanned by H along metric defined by R = chol(r_cov)
# returns this projection matrix multiplied by r_cov, which is symmetric
oblique_proj = function(H, R) {
    H_tilde = backsolve(R, H, transpose = TRUE)
    tcrossprod(H %*% solve(qr.R(qr(H_tilde))))
}
