#' Estimate ecological quantities
#'
#' Produces estimates of overall conditional means from a fitted
#' ecological inference model or Riesz representer.
#' If both a regression model and a Riesz representer are provided, a debiased
#' machine learning (DML) estimate is produced.
#'
#' @param regr A fitted regression model, from [ei_ridge()], or another kind
#'    of regression model wrapped with [ei_wrap_model()].
#'    If `riesz` is not provided and `regr` is an [ei_riesz()] object, then
#'    `riesz` will be set to the value of `regr` and `regr` will be set to
#'    `NULL`. This is so users can call this function as
#'    `ei_est(<riesz>, data = <data>)`.
#' @param riesz A fitted Riesz representer, from [ei_riesz()], or a matrix of
#'    Riesz weights
#' @param data The data frame, matrix, or [ei_spec] object that was used to fit
#'   the regression or Riesz representer.
#' @param total <[`tidy-select`][dplyr::select]> A variable containing the total
#'   number of observations in each aggregate unit. For example, the column
#'   containing the total number of voters. Required if `data` is not an
#'   [ei_spec()] object and `riesz` is not provided.
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
#'                total = pres_total, covariates = c(pop_urban, farm))
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
ei_est = function(regr=NULL, riesz=NULL, data, total, outcome=NULL, conf_level=FALSE) {
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

    if (missing(total)) {
        if (inherits(riesz, "ei_riesz")) {
            total = riesz$blueprint$ei_n
        }
        if (inherits(data, "ei_spec")) {
            total = attr(data, "ei_n")
        }
    }
    w = check_make_weights(!!enquo(total), data, n)
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
    # TODO check
    # riesz = riesz * weights
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
    if (inherits(regr, "ei_wrapped")) {
        return(regr)
    }
    if (!inherits(regr, c("ei_ridge", "ei_model"))) {
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

    list(yhat=regr$fitted, preds=preds, x=x)
}

#' Wrap another predictive model for use in `ei_est`
#'
#' Stores additional data and attributes on a generic model class so that it
#' can be used as the `regr` argument to [ei_est()]. Given the wide variety of
#' model classes, there is no guarantee this function will work. However, most
#' model classes supporting a [fitted()] and [predict()] method will work as long
#' as there is no transformation of the predictor variables as part of the model
#' formula or fitting.
#'
#' @param x A model object, supporting [fitted()] and [predict()] generics.
#' @param data A data frame or matrix containing the data used to fit the model,
#'   or an [ei_spec()] object (recommended). If the latter, then the
#'   `predictors` and `outcome` arguments are ignored and need not be provided.
#' @inheritParams ei_spec
#'
#' @returns An `ei_wrapped` object, which has the information required to use
#'   the provided `x` with [ei_est()].
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal, pres_total,
#'                covariates = c(pop_urban, farm))
#'
#' m = lm(pres_ind_wal ~ 0 + white + black + other + pop_urban + farm, spec)
#' m_wrap = ei_wrap_model(m, spec)
#' print(m_wrap)
#'
#' ei_est(m_wrap, data = spec)
#'
#' @export
ei_wrap_model <- function(x, data, predictors = NULL, outcome = NULL) {
    if (inherits(data, "ei_spec")) {
        predictors = attr(data, "ei_x")
        outcome = attr(data, "ei_y")
    } else {
        predictors = try_fetch(
            names(eval_select(enquo(predictors), data, allow_empty=FALSE)),
            error = function(cnd) cli_abort("Predictor specification failed.", parent=cnd)
        )
        outcome =  try_fetch(
            names(eval_select(enquo(outcome), data, allow_empty=FALSE)),
            error = function(cnd) cli_abort("Outcome specification failed.", parent=cnd)
        )
    }

    # find predict() and its newdata argument
    fn_pred = NULL
    for (cl in class(x)) {
        fn_pred = utils::getS3method("predict", cl, optional=TRUE)
        if (!is.null(fn_pred))
            break
    }
    if (is.null(fn_pred)) {
        cli_abort("No {.fn predict} method for an object
                  of class {.cls {class(x)}}", call = parent.frame())
    }

    args = rlang::fn_fmls_names(fn_pred)
    arg = args[grepl("new.?data", args)]
    if (length(arg) != 1) {
        cli_abort("No {.arg new_data}, {.arg newdata}, or similar argument to the
                  {.fn predict} method for an object of class {.cls {class(x)}}",
                  call = parent.frame())
    }

    preds = list()
    for (group in predictors) {
        data_copy = data
        data_copy[, predictors] = 0
        data_copy[, group] = 1

        argl = list(x, data_copy)
        names(argl) = c("", arg)
        preds[[group]] = as.matrix(do.call(fn_pred, argl))
    }

    out = list(
        y = as.matrix(data[, outcome]),
        yhat = fitted(x),
        preds = preds,
        x = as.matrix(data[, predictors]),
        blueprint = list(ei_x = predictors),
        classes = class(x)
    )
    class(out) = "ei_wrapped"

    if (is.null(out$yhat)) {
        cli_abort("Unable to extract fitted values for {.cls class(x)} model.")
    }

    out
}


# Types --------

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

#' @export
print.ei_wrapped <- function(x, ...) {
    cat_line(format_inline("A wrapped {.cls {x$classes}} model with {length(x$yhat)} observations"))
}
