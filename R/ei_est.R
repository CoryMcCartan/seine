#' Estimate ecological quantities
#'
#' Produces estimates of overall conditional means from a fitted
#' ecological inference model or Riesz representer.
#' If both a regression model and a Riesz representer are provided, a debiased
#' machine learning (DML) estimate is produced.
#'
#' @param regr A fitted regression model, from [ei_ridge()].
#'    If `riesz` is not provided and `regr` is an [ei_riesz()] object, then
#'    `riesz` will be set to the value of `regr` and `regr` will be set to
#'    `NULL`. This is so users can call this function as
#'    `ei_est(<riesz>, data = <data>)`.
#' @param riesz A fitted Riesz representer, from [ei_riesz()], or a matrix of
#'    Riesz weights
#' @param data The data frame, matrix, or [ei_spec] object that was used to fit
#'   the regression or Riesz representer.
#' @param weights <[`data-masking`][rlang::args_data_masking]> A vector of unit
#'   weights. In general these should be the count of the total number of
#'   individuals in each unit. Required by default. To force constant weights,
#'   you can provide `weights=FALSE`.
#' @param outcome <[`data-masking`][rlang::args_data_masking]> A vector or
#'   matrix of outcome variables. Only required if both `riesz` is provided
#'   alone (without `regr`) and `data` is not an [ei_spec] object.
#' @param conf_level A numeric specifying the level for confidence intervals.
#'   If `FALSE` (the default), no confidence intervals are calculated. Standard
#'   errors are always returned.
#'
#' @returns A data frame with estimates.
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
#'                weights = pres_total, covariates = c(pop_urban, farm))
#'
#' m = ei_ridge(spec)
#' rr = ei_riesz(spec, penalty = m$penalty)
#'
#' ei_est(regr = m, data = spec, conf_level = 0.95) # Plug-in estimate
#' ei_est(riesz = rr, data = spec) # Weighted (Riesz) estimate
#' est = ei_est(regr = m, riesz = rr, data = spec) # Double/debiased ML estimate
#' as.matrix(est)
#' as.matrix(est, se=TRUE)
#' vcov(est)[1:4, 1:4]
#'
#' @export
ei_est = function(regr=NULL, riesz=NULL, data, weights, outcome=NULL, conf_level=FALSE) {
    if (is.null(regr) && is.null(riesz)) {
        cli_abort("At least one of {.arg regr} or {.arg riesz} must be provided.")
    }
    if (is.null(riesz) && inherits(regr, "ei_riesz")) {
        riesz = regr
        regr = NULL
    }

    y = est_check_outcome(regr, data, !!enquo(outcome))
    n = nrow(y)
    n_y = ncol(y)

    if (missing(weights) && inherits(data, "ei_spec")) {
        weights = attr(data, "ei_wgt")
    }
    w = check_make_weights(!!enquo(weights), data, n)
    w = w / mean(w)

    # build predictions and RR
    riesz = est_check_riesz(riesz, data, w, n, regr)
    rl = est_check_regr(regr, data, n, colnames(riesz), n_y)

    n_x = ncol(riesz)
    xc = colnames(riesz)
    eif = matrix(nrow=n, ncol=n_y*n_x) # x varies faster than y
    for (i in seq_len(n_x)) {
        plugin = rl$preds[[xc[i]]] * rl$x[, i] * w / mean(rl$x[, i] * w)
        wtd = riesz[, i] * (y - rl$yhat)
        eif[, i + (seq_len(n_y) - 1)*n_x] = plugin + wtd
    }
    est = colMeans(eif)
    vcov = crossprod(shift_cols(eif, est)) / (n - 1)^2
    se = sqrt(diag(vcov))

    out = tibble::new_tibble(list(
        predictor = rep(xc, n_y),
        outcome = rep(colnames(y), each=n_x),
        estimate = est,
        std.error = se
    ), class="ei_est")
    if (!isFALSE(conf_level)) {
        alpha = (1 - conf_level) / 2
        out$conf.low = out$estimate + qnorm(alpha) * out$std.error
        out$conf.high = out$estimate - qnorm(alpha) * out$std.error
    }

    rownames(vcov) = colnames(vcov) = c(outer(xc, colnames(y), paste, sep=":"))
    attr(out, "vcov") = vcov

    out
}

est_check_outcome = function(regr, data, outcome) {
    if (!is.null(regr)) {
        y = regr$y
    } else if (inherits(data, "ei_spec")) {
        y = as.matrix(data[attr(data, "ei_y")])
    } else {
        quo = enquo(outcome)
        y = as.matrix(eval_tidy(!!quo, data))
        if (ncol(y) == 1) {
            colnames(y) = rlang::quo_label(quo)
        }
    }
    if (is.null(y)) {
        cli_abort("The {.arg outcome} argument is required when {.arg riesz}
                   is provided without {.arg regr} and {.arg data} is not
                   an {.cls ei_spec} object.", call=parent.frame())
    }

    y
}

est_check_riesz = function(riesz, data, weights, n, regr) {
    if (is.null(riesz)) {
        xcols = regr$blueprint$ei_x # will be NULL if regr is wrong type
        riesz = matrix(1/n, nrow=n, ncol=length(xcols))
        colnames(riesz) = xcols
    } else if (inherits(riesz, "ei_riesz")) {
        riesz = riesz$weights
    } else if (!is.matrix(riesz)) {
        cli_abort("{.arg riesz} must be an {.cls ei_ridge} object or a matrix.",
                  call=parent.frame())
    }

    if (nrow(riesz) != n) {
        cli_abort("The number of weights in {.arg riesz} ({nrow(riesz)}) must
                  match the number of observations ({n}).", call=parent.frame())
    }

    riesz = riesz * weights # TODO FIX
    riesz = scale_cols(riesz, 1 / colMeans(riesz))
    riesz
}

est_check_regr = function(regr, data, n, xcols, n_y) {
    if (is.null(regr)) {
        preds = lapply(xcols, function(.) matrix(0, nrow=n, ncol=n_y))
        names(preds) = xcols
        x = matrix(1, nrow=n, ncol=length(xcols))
        return(list(yhat=preds[[1]], preds=preds, x=x, z=NULL))
    }
    if (!inherits(regr, "ei_ridge") && !inherits(regr, "ei_model")) {
        cli_abort("{.arg regr} must be a {.cls ei_ridge} object.",
                  call=parent.frame())
    }

    data = hardhat::forge(data, regr$blueprint)$predictors
    idx_x = match(regr$blueprint$ei_x, colnames(data))
    z = as.matrix(data[, -idx_x, drop=FALSE])
    x = data[, idx_x, drop=FALSE]
    n_x = length(xcols)
    p = ncol(z)

    # normalize and check
    z = shift_cols(z, regr$z_shift)
    z = scale_cols(z, regr$z_scale)
    if (any(is.na(z)))
        cli_abort("Missing values found in covariates.", call=parent.frame())

    preds = list()
    for (group in seq_along(xcols)) {
        use = n_x + p*(group-1) + seq_len(p)
        preds[[xcols[group]]] = shift_cols(z %*% regr$coef[use, ],
                                           regr$coef[group, ] * -regr$int_scale)
    }

    list(yhat=regr$fitted, preds=preds, x=x, z=z)
}


# Type --------

#' @describeIn ei_est Format estimates or standard errors as a matrix
#' @param x,object An object of class `ei_est`
#' @param se If `TRUE`, return standard errors instead of estimates.
#' @param ... Additional arguments (ignored)
#' @export
as.matrix.ei_est <- function(x, se=FALSE, ...) {
    y = if (isFALSE(se)) x$estimate else x$std.error
    xlab = unique(x$predictor)
    ylab = unique(x$outcome)
    matrix(y, nrow=length(xlab), ncol=length(ylab), dimnames=list(xlab, ylab))
}

#' @describeIn ei_est Extract full covariance matrix of estimates
#' @export
vcov.ei_est <- function(object, ...) {
    attr(object, "vcov")
}
