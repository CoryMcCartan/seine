#' Target an ecological regression to respect constraints on the outcome
#'
#' Updates a fitted ecological regression using Targeted Minimum Loss-based
#' Estimation (TMLE). The returned targeted regression can be supplied to
#' [ei_est()] together with the same Riesz representer to produce a targeted
#' estimate that respects the provided outcome bounds and optional sum-to-one
#' constraint.
#'
#' @param regr A fitted regression model, from [ei_ridge()], or another kind
#'   of regression model wrapped with [ei_wrap_model()].
#' @param riesz A fitted Riesz representer from [ei_riesz()].
#' @param data The data frame, matrix, or [ei_spec] object used to fit the
#'   regression and Riesz representer.
#' @param bounds A vector `c(min, max)` of bounds for the outcome. By default,
#'   uses the bounds stored in `regr`, if available; otherwise infers bounds
#'   from the outcome variable. Set to `FALSE` for unbounded targeting.
#' @param sum_one If `TRUE`, constrain outcomes to sum to one. By default, uses
#'   the value stored in `regr`, if available; otherwise infers the constraint
#'   from the outcome variable and bounds.
#'
#' @returns An `ei_targeted` object, which is also an `ei_wrapped` object and
#'   can be supplied as the `regr` argument to [ei_est()].
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
#'                total = pres_total, covariates = c(pop_urban, farm))
#' m = ei_ridge(spec)
#' rr = ei_riesz(spec, penalty = m$penalty)
#'
#' m_targeted = ei_target(m, rr, spec, bounds = 0:1, sum_one = TRUE)
#' ei_est(m_targeted, rr, spec)
#' @export
ei_target = function(regr, riesz=NULL, data, bounds=NULL, sum_one=NULL) {
    if (!inherits(riesz, "ei_riesz")) {
        cli_abort("{.arg riesz} must be an {.cls ei_riesz} object.",
                  call=parent.frame())
    }

    y = est_check_outcome(regr, data, NULL)
    n = nrow(y)
    n_y = ncol(y)
    rm = est_check_riesz(riesz, data, rep(1, n), n, regr)
    rl = est_check_regr(regr, data, n, colnames(rm), n_y, vcov=FALSE)
    pred_vcov_u = NULL
    if (inherits(regr, "ei_ridge") && !is.null(regr$vcov_u)) {
        pred_vcov_u = est_check_regr(regr, data, n, colnames(rm), n_y, vcov=TRUE)$vcov_u
    }

    if (is.null(sum_one) && !is.null(regr$blueprint$sum_one) && !isFALSE(bounds)) {
        sum_one = regr$blueprint$sum_one
    }
    if (is.null(bounds) && !is.null(regr$blueprint$bounds)) {
        bounds = regr$blueprint$bounds
    }
    bounds = check_bounds(bounds, y, clamp=1e-8)
    if (is.null(sum_one) && all(bounds == c(0, 1))) {
        sum_one = isTRUE(all.equal(rowSums(y), rep(1, nrow(y))))
    }

    rm_cf = tmle_riesz_cf(riesz, data)
    rl = tmle_target(y, rl, rm, rm_cf, bounds, sum_one)

    w = regr$blueprint$ei_wgt
    if (is.null(w)) {
        w = rep(1, n)
    }

    structure(list(
        y = y,
        yhat = rl$yhat,
        sigma2 = colMeans((y - rl$yhat)^2 * w / mean(w)),
        eps = rl$eps,
        vcov_u = regr$vcov_u,
        pred_vcov_u = pred_vcov_u,
        r2 = diag(as.matrix(cor(rl$yhat, y)^2)),
        preds = rl$preds,
        x = rl$x,
        blueprint = list(
            ei_x = names(rl$preds),
            ei_wgt = w,
            bounds = bounds,
            sum_one = sum_one
        ),
        classes = if (inherits(regr, "ei_wrapped")) regr$classes else class(regr)
    ), class=c("ei_targeted", "ei_wrapped"))
}

tmle_target = function(y, rl, rm, rm_cf, bounds, sum_one) {
    fluctuation = tmle_fluctuation(bounds, sum_one)
    rl$yhat = fluctuation$project(rl$yhat)
    rl$preds = lapply(rl$preds, fluctuation$project)

    eps = fluctuation$fit_eps(y, rl$yhat, rm)
    dimnames(eps) = list(colnames(rm), colnames(y))
    rl$yhat = fluctuation$update(rl$yhat, rm, eps)

    score = crossprod(rm, y - rl$yhat) / nrow(y)
    if (max(abs(score)) > 1e-3) {
        cli_warn("Targeting step may not have converged properly.")
    }

    for (group in names(rl$preds)) {
        rl$preds[[group]] = fluctuation$update(rl$preds[[group]], rm_cf[[group]], eps)
    }
    rl$eps = eps
    rl
}

# Construct the counterfactual 'clever covariates' for a fitted Riesz object.
tmle_riesz_cf = function(riesz, data) {
    data = hardhat::forge(data, riesz$blueprint)$predictors
    xcols = riesz$blueprint$ei_x
    idx_x = match(xcols, colnames(data))
    z = as.matrix(data[, -idx_x, drop = FALSE])
    p = ncol(z)
    n_x = length(xcols)

    z = shift_cols(z, riesz$z_shift)
    z = scale_cols(z, riesz$z_scale)
    if (anyNA(z)) {
        cli_abort("Missing values found in covariates.", call=parent.frame())
    }
    z = cbind(riesz$int_scale, z)
    w = riesz$blueprint$ei_wgt / mean(riesz$blueprint$ei_wgt)
    w = w %o% (riesz$int_scale / colMeans(riesz$weights))

    rm_cf = list()
    for (group in seq_along(xcols)) {
        use = c(group, n_x + p*(group - 1) + seq_len(p))
        rm_cf[[xcols[group]]] = w * (z %*% riesz$coef[use, ])
    }
    rm_cf
}



# Select a fluctuation submodel which preserves the requested outcome support.
tmle_fluctuation = function(bounds, sum_one) {
    if (isTRUE(sum_one)) {
        if (!isTRUE(all.equal(bounds, c(0, 1)))) {
            cli_abort(
                "TMLE with {.arg sum_one = TRUE} currently supports only {.code bounds = c(0, 1)}.",
                call = parent.frame()
            )
        }
        list(
            project = tmle_project_simplex,
            fit_eps = function(y, q, H) tmle_fit_eps(y, q, H, tmle_loss_multinomial),
            update = tmle_update_softmax
        )
    } else if (all(is.finite(bounds))) {
        scale = bounds[2] - bounds[1]
        list(
            project = function(q) tmle_project_interval(q, bounds[1], bounds[2]),
            fit_eps = function(y, q, H) {
                tmle_fit_eps(
                    (y - bounds[1]) / scale,
                    (q - bounds[1]) / scale,
                    H,
                    tmle_loss_binomial
                )
            },
            update = function(q, H, eps) tmle_update_logistic(q, H, eps, bounds)
        )
    } else if (is.finite(bounds[1])) {
        list(
            project = function(q) tmle_project_lower(q, bounds[1]),
            fit_eps = function(y, q, H) {
                tmle_fit_eps(y - bounds[1], q - bounds[1], H, tmle_loss_poisson)
            },
            update = function(q, H, eps) tmle_update_log(q, H, eps, bounds[1])
        )
    } else if (is.finite(bounds[2])) {
        list(
            project = function(q) tmle_project_upper(q, bounds[2]),
            fit_eps = function(y, q, H) {
                tmle_fit_eps(bounds[2] - y, bounds[2] - q, H, tmle_loss_poisson)
            },
            update = function(q, H, eps) tmle_update_log(q, H, eps, bounds[2], upper = TRUE)
        )
    } else {
        if (isTRUE(sum_one)) {
            cli_abort("TMLE cannot enforce {.arg sum_one = TRUE} for an unbounded outcome.",
                      call = parent.frame())
        }
        list(
            project = identity,
            fit_eps = tmle_fit_eps_gaussian,
            update = function(q, H, eps) q + H %*% eps
        )
    }
}

# Fit an unbounded least-squares fluctuation after column-normalizing H.
tmle_fit_eps_gaussian = function(y, q, H) {
    H_scale = pmax(sqrt(colMeans(H^2)), 1e-8)
    H_norm = scale_cols(H, 1 / H_scale)
    udv = svd(H_norm)
    tol = max(dim(H_norm)) * max(udv$d) * .Machine$double.eps
    d_inv = ifelse(udv$d > tol, 1 / udv$d, 0)
    eps = udv$v %*% (d_inv * crossprod(udv$u, y - q))
    sweep(eps, 1, H_scale, "/")
}


# Fit a J by K fluctuation parameter matrix after column-normalizing H.
tmle_fit_eps = function(y, q, H, loss) {
    n_x = ncol(H)
    n_y = ncol(y)
    H_scale = pmax(sqrt(colMeans(H^2)), 1e-8)
    H_norm = scale_cols(H, 1 / H_scale)

    nllik = function(eps_vec) {
        eps = matrix(eps_vec, n_x, n_y)
        loss(y, q, H_norm, eps)
    }
    opt = optim(
        par = rep(0, n_x * n_y),
        fn = nllik,
        method = "L-BFGS-B",
        control = list(maxit = 1000, factr = 1e3)
    )
    matrix(opt$par, n_x, n_y) / H_scale
}

tmle_loss_multinomial = function(y, q, H, eps) {
    linpred = log(q) + H %*% eps
    max_lp = linpred[cbind(seq_len(nrow(linpred)), max.col(linpred))]
    linpred = linpred - max_lp
    sum(log(rowSums(exp(linpred))) - rowSums(y * linpred))
}

tmle_loss_binomial = function(y, q, H, eps) {
    linpred = qlogis(q) + H %*% eps
    sum(pmax(linpred, 0) + log1p(exp(-abs(linpred))) - y * linpred)
}

tmle_loss_poisson = function(y, q, H, eps) {
    linpred = log(q) + H %*% eps
    sum(exp(linpred) - y * linpred)
}

tmle_update_softmax = function(q, H, eps) {
    linpred = log(q) + H %*% eps
    max_lp = linpred[cbind(seq_len(nrow(linpred)), max.col(linpred))]
    preds = exp(linpred - max_lp)
    preds / rowSums(preds)
}

tmle_update_logistic = function(q, H, eps, bounds) {
    scale = bounds[2] - bounds[1]
    bounds[1] + scale * plogis(qlogis((q - bounds[1]) / scale) + H %*% eps)
}

tmle_update_log = function(q, H, eps, bound, upper = FALSE) {
    q = if (upper) bound - q else q - bound
    q = q * exp(H %*% eps)
    if (upper) bound - q else bound + q
}

# Duchi et al. (2008) fast simplex projection onto its delta-interior.
tmle_project_simplex = function(x, delta = 1e-8) {
    U = t(apply(x - delta, 1, function(z) sort.int(z, decreasing = TRUE, method = "quick")))
    z = 1 - ncol(x) * delta

    L = diag(ncol(x))
    L[upper.tri(L)] = 1
    cssv = U %*% L
    idxs = rep(seq_len(ncol(x)), each = nrow(x))
    rho = rowSums(U - (cssv - z) / idxs > 0)

    cssv_rho = cssv[cbind(seq_len(nrow(x)), rho)]
    theta = (cssv_rho - z) / rho
    x = x - theta
    x[x < 0] = 0
    x + delta
}

tmle_project_interval = function(q, lower, upper, delta = 1e-8) {
    delta = min(delta, (upper - lower) / 4)
    pmin(pmax(q, lower + delta), upper - delta)
}

tmle_project_lower = function(q, lower, delta = 1e-8) {
    pmax(q, lower + delta)
}

tmle_project_upper = function(q, upper, delta = 1e-8) {
    pmin(q, upper - delta)
}

#' @export
print.ei_targeted = function(x, ...) {
    cat_line(format_inline(paste0(
        "A targeted {.cls {x$classes[1]}} model with {ncol(x$y)} outcome{?s}, ",
        "{length(x$blueprint$ei_x)} group{?s}, and {nrow(x$yhat)} observations"
    )))
    bounds = x$blueprint$bounds
    if (any(is.finite(bounds))) {
        sumt1 = if (isTRUE(x$blueprint$sum_one)) " and constrained to sum to 1" else ""
        pl = if (ncol(x$y) > 1) "s" else ""
        cat_line("With outcome", pl, " bounded in (", bounds[1], ", ", bounds[2], ")", sumt1)
    }
}

#' @export
fitted.ei_targeted = function(object, ...) {
    object$yhat
}

#' @export
residuals.ei_targeted = function(object, ...) {
    object$y - object$yhat
}

#' @export
vcov.ei_targeted = function(object, ...) {
    object$vcov_u
}

#' @export
summary.ei_targeted = function(object, ...) {
    cat_line("Fitted values:")
    print(summary((object$yhat)))
    cat_line()
    cat_line("R-squared by outcome:")
    print(object$r2)
}

#' @export
weights.ei_targeted = function(object, normalize = TRUE, ...) {
    w = object$blueprint$ei_wgt
    if (is.null(w)) w = rep(1, nrow(object$y))
    if (isTRUE(normalize)) w / mean(w) else w
}
