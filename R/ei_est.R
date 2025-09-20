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
#' @param total <[`tidy-select`][tidyselect::select_helpers]> A variable
#'   containing the total number of observations in each aggregate unit. For
#'   example, the column containing the total number of voters. Required if
#'   `data` is not an [ei_spec()] object and `riesz` is not provided.
#' @param subset <[`data-masking`][rlang::args_data_masking]> An optional
#'   indexing vector describing the subset of units over which to calculate
#'   estimates.
#' @param outcome <[`data-masking`][rlang::args_data_masking]> A vector or
#'   matrix of outcome variables. Only required if both `riesz` is provided
#'   alone (without `regr`) and `data` is not an [ei_spec] object.
#' @param conf_level A numeric specifying the level for confidence intervals.
#'   If `FALSE` (the default), no confidence intervals are calculated. Standard
#'   errors are always returned.
#' @param use_student If `TRUE`, use construct confidence intervals from a
#'   Student-_t_ distribution, which may improve coverage properties in
#'   small samples.
#'
#' @returns A data frame with estimates. It has class `ei_est`, supporting
#'   several methods, and two additional attributes: `vcov`, containing the
#'   estimated covariance matrix for the estimates, and `n`, containing the
#'   number of aggregate units used in estimation (the number of rows in
#'   `data`).
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
#'                total = pres_total, covariates = c(state, pop_urban, farm))
#'
#' m = ei_ridge(spec)
#' rr = ei_riesz(spec, penalty = m$penalty)
#'
#' ei_est(regr = m, data = spec, conf_level = 0.95) # Plug-in estimate
#' ei_est(riesz = rr, data = spec) # Weighted (Riesz) estimate
#' est = ei_est(regr = m, riesz = rr, data = spec) # Double/debiased ML estimate
#' as.matrix(est)
#' as.matrix(est, "std.error")
#' vcov(est)[1:4, 1:4]
#'
#' est = ei_est(m, rr, data = spec, subset = (state == "Alabama"))
#' as.matrix(est)
#' nobs(est)
#' @export
ei_est = function(regr=NULL, riesz=NULL, data, total, subset=NULL,
                  outcome=NULL, conf_level=FALSE, use_student=TRUE) {
    if (is.null(regr) && is.null(riesz)) {
        cli_abort("At least one of {.arg regr} or {.arg riesz} must be provided.")
    }
    if (is.null(riesz) && inherits(regr, "ei_riesz")) {
        riesz = regr
        regr = NULL
    }

    y = est_check_outcome(regr, data, enquo(outcome))
    n = nrow(y)
    n_y = ncol(y)

    subset = eval_tidy(enquo(subset), data)
    if (!is.null(subset)) {
        if (is.logical(subset)) {
            subset = 1*subset
        } else if (is.numeric(subset)) {
            subset = 1*(seq_len(n) == subset)
        } else {
            cli_abort("The {.arg subset} argument must be a logical or integer
                       vector.", call=parent.frame())
        }
        if (any(is.na(subset))) {
            cli_abort("The {.arg subset} argument must not contain missing values.",
                      call=parent.frame())
        }
    } else {
        subset = rep(1, n)
    }

    if (missing(total)) {
        if (inherits(riesz, "ei_riesz")) {
            total = riesz$blueprint$ei_n
        }
        if (inherits(data, "ei_spec")) {
            total = attr(data, "ei_n")
        }
    }
    w = check_make_weights(!!enquo(total), data, n)
    w = subset * w / mean(subset * w)

    # save data for sensitivity analysis
    if (inherits(riesz, "ei_riesz") && inherits(regr, "ei_ridge")) {
        sens_s = sqrt(riesz$nu2 %o% regr$sigma2)
    } else {
        sens_s = NULL
    }

    # build predictions and RR
    riesz = est_check_riesz(riesz, data, w, n, regr)
    rl = est_check_regr(regr, data, n, colnames(riesz), n_y)

    n_x = ncol(riesz)
    xc = names(rl$preds)
    eif = matrix(nrow=n, ncol=n_y*n_x) # x varies faster than y
    for (i in seq_len(n_x)) {
        plugin = rl$preds[[xc[i]]] * rl$x[, i] * w / mean(rl$x[, i] * w)
        wtd = riesz[, i] * (y - rl$yhat)
        eif[, i + (seq_len(n_y) - 1)*n_x] = plugin + wtd
    }
    est = colMeans(eif)
    vcov = crossprod(shift_cols(eif, est)) / (sum(subset) - 1)^2
    se = sqrt(diag(vcov))

    out = tibble::new_tibble(list(
        predictor = rep(xc, n_y),
        outcome = rep(colnames(y), each=n_x),
        estimate = est,
        std.error = se
    ), class="ei_est")
    if (!isFALSE(conf_level)) {
        alpha = (1 - conf_level) / 2
        crit = if (isFALSE(use_student)) {
            qnorm(alpha)
        } else {
            qt(alpha, df = n - 1)
        }
        out$conf.low = out$estimate + crit * out$std.error
        out$conf.high = out$estimate - crit * out$std.error
    }

    rownames(vcov) = colnames(vcov) = c(outer(xc, colnames(y), paste, sep=":"))
    attr(out, "vcov") = vcov
    attr(out, "n") = sum(subset)
    attr(out, "sens_s") = sens_s
    attr(out, "bounds_inf") = ei_bounds(NULL, rl$y)

    out
}

est_check_outcome = function(regr, data, quo_outcome) {
    if (!is.null(regr)) {
        y = regr$y
    } else if (inherits(data, "ei_spec")) {
        y = as.matrix(data[attr(data, "ei_y")])
    } else {
        y = eval_tidy(quo_outcome, data)
        if (!is.null(y)) {
            y = as.matrix(y)
            if (ncol(y) == 1) {
                colnames(y) = rlang::as_label(quo_outcome)
            }
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
        if (length(xcols) == 1) {
            xcols = c(xcols, ".other")
        }
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
    riesz = scale_cols(riesz, 1 / colMeans(riesz))
    riesz
}

est_check_regr = function(regr, data, n, xcols, n_y, sd = FALSE) {
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
    xcols = regr$blueprint$ei_x
    idx_x = match(xcols, colnames(data))
    z = as.matrix(data[, -idx_x, drop=FALSE])
    x = pull_x(data, idx_x)
    xcols = colnames(x)
    n_x = length(xcols)
    p = ncol(z)

    # normalize and check
    z = shift_cols(z, regr$z_shift)
    z = scale_cols(z, regr$z_scale)
    if (any(is.na(z)))
        cli_abort("Missing values found in covariates.", call=parent.frame())
    z = cbind(regr$int_scale, z)

    preds = list()
    sds = if (sd) matrix(nrow = n, ncol = n_x^2) else
    for (group in seq_along(xcols)) {
        use = c(group, n_x + p*(group-1) + seq_len(p))
        preds[[xcols[group]]] = z %*% regr$coef[use, ]

        if (sd) {
            for (g2 in seq_len(group)) {
                use2 = c(g2, n_x + p*(g2-1) + seq_len(p))
                sds[, (group - 1)*n_x + g2] = rowSums(z * (z %*% regr$vcov_u[use, use2]))
                sds[, (g2 - 1)*n_x + group] = sds[, (group - 1)*n_x + g2]
            }
        }
    }

    list(yhat=regr$fitted, preds=preds, sds=sds, x=x)
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
#' @param ... Additional arguments passed to the [predict()] method.
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
#' # Note: this is not a model recommended for valid ecological inference!
#' m = suppressWarnings(
#'     glm(pres_ind_wal ~ 0 + vap_white + vap_black + vap_other + pop_urban + farm,
#'         data = spec, family = "binomial")
#' )
#' m_wrap = ei_wrap_model(m, spec, type = "response")
#' print(m_wrap)
#'
#' ei_est(m_wrap, data = spec) # notice all estimates nonnegative
#'
#' @export
ei_wrap_model <- function(x, data, predictors = NULL, outcome = NULL, ...) {
    if (inherits(x, "ei")) {
        return(wrap_king_ei(x))
    }
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

        argl = list(x, data_copy, ...)
        names(argl) = c("", arg, names(argl)[-1:-2])
        preds[[group]] = as.matrix(do.call(fn_pred, argl))
    }

    y = as.matrix(data[, outcome])
    yhat = as.matrix(fitted(x))
    sigma2 = colMeans((y - yhat)^2)
    out = list(
        y = y,
        yhat = yhat,
        sigma2 = sigma2,
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

wrap_king_ei <- function(obj) {
    y = as.matrix(obj$t)
    n = length(y)
    x = cbind(.x = obj$x, .other = 1 - obj$x)

    if (!isTRUE(all.equal(obj$Zb, obj$Zw))) {
        cli_abort("Covariates for each group, {.arg Zb} and {.arg Zw}, must match.",
                   call = parent.frame())
    }

    # calculate CEF manually
    # $betaw and $betaw combine the CEF with the residuals
    # this approach also obviates any sampling, by using `ep_moments`
    sds = exp(obj$phi[3:4])
    Sigma = (diag(2)*(1 - tanh(obj$phi[5])) + tanh(obj$phi[5]))
    Sigma = diag(sds) %*% Sigma %*% diag(sds)
    L = t(chol(Sigma))
    b_loc = obj$phi[1:2] * (sds^2 + 0.25) + 0.5

    z = shift_cols(obj$Zb, colMeans(obj$Zb))
    eta_raw = shift_cols(z %*% matrix(obj$phi[-1:-5], ncol=2), -b_loc)
    eta = matrix(nrow=n, ncol=2)
    for (i in seq_len(n)) {
        eta[i, ] = R_ep_moments(eta_raw[i, ], L, numeric(0), 0, 1e-7)[[2]]
    }

    yhat = rowSums(eta * x)

    structure(list(
        y = y,
        yhat = yhat,
        sigma2 = var(y - yhat),
        preds = list(.x = eta[, 1, drop=FALSE], .other = eta[, 2, drop=FALSE]),
        x = x,
        blueprint = list(ei_x = colnames(x)),
        classes = "ei"
    ), class = "ei_wrapped")
}


# Types --------

#' @describeIn ei_est Format estimates, standard errors, or other columns as a matrix.
#' @param x,object An object of class `ei_est`
#' @param which Which column of `ei_est` to convert to a matrix. For example,
#'   pass `which="std.error"` to return standard errors instead of estimates.
#'   Partial matching supported.
#' @param ... Additional arguments (ignored)
#' @export
as.matrix.ei_est <- function(x, which="estimate", ...) {
    nms = setdiff(names(x), c("predictor", "outcome"))
    col = nms[pmatch(which, nms)]
    if (is.na(col)) {
        cli_abort(c("Unknown column {.val {which}} in {.cls ei_est} object.",
                    ">"="Available columns: {.val {nms}}"))
    }
    y = x[[col]]

    xlab = unique(x$predictor)
    ylab = unique(x$outcome)
    out = matrix(nrow=length(xlab), ncol=length(ylab), dimnames=list(xlab, ylab))
    out[cbind(x$predictor, x$outcome)] = y
    names(dimnames(out)) = c("predictor", "outcome")
    out
}

#' @describeIn ei_est Extract full covariance matrix of estimates
#' @export
vcov.ei_est <- function(object, ...) {
    attr(object, "vcov")
}

#' @describeIn ei_est Extract number of units covered by estimates
#' @export
nobs.ei_est <- function(object, ...) {
    attr(object, "n")
}

#' @export
print.ei_wrapped <- function(x, ...) {
    cat_line(format_inline("A wrapped {.cls {x$classes}} model with {length(x$yhat)} observations"))
}
