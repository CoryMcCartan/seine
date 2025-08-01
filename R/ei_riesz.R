# User-facing model functions --------------------------------------------------

#' Estimate Riesz representer for ecological inference
#'
#' Fits a penalized Riesz regression for ecological inference, allowing for
#' overall estimates of conditional means using [ei_est()].
#'
#' The regression is calculated using the singular value decomposition.
#'
#' @param penalty The ridge penalty (a non-negative scalar), which must be
#'   specified. Recommended value is the same penalty used in [ei_ridge()],
#'   which is stored in the `penalty` entry of the fitted model object.
#' @inheritSection ei_ridge Weights
#'
#' @inheritParams ei_ridge
#' @inheritParams ei_spec
#'
#' @returns An `ei_riesz` object.
#'
#' @examples
#' data(elec_1968)
#'
#' # Recommended: get ridge penalty from ei_ridge()
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
#'                total = pres_total, covariates = c(pop_urban, farm))
#' m = ei_ridge(spec)
#'
#' ei_riesz(spec, penalty = m$penalty)
#'
#' rr = ei_riesz(~ vap_white + vap_black + vap_other | pop_urban + farm,
#'               data = elec_1968, total = pres_total, penalty = m$penalty)
#' summary(rr)
#'
#' # Examine the weights and check they have mean 1
#' head(weights(rr, group = "vap_black"))
#' colMeans(weights(rr))
#'
#' mean(elec_1968$pres_ind_wal * weights(rr, "vap_white"))
#' @export
ei_riesz <- function(x, ...) {
    UseMethod("ei_riesz")
}


#' @export
#' @rdname ei_riesz
ei_riesz.formula <- function(formula, data, total, weights, penalty, scale=TRUE, ...) {
    f_lhs(formula) = NULL
    forms = ei_forms(formula)
    form_preds = terms(rlang::new_formula(lhs=NULL, rhs=forms$predictors))
    form_combined = rlang::new_formula(forms$outcome, expr(!!forms$predictors + !!forms$covariates))

    bp = hardhat::new_default_formula_blueprint(
        intercept = FALSE,
        composition = "matrix",
        indicators = "one_hot",
        ei_x = attr(form_preds, "term.labels"),
        ei_n = check_make_weights(!!enquo(total), data),
        ei_wgt = check_make_weights(!!enquo(weights), data, arg="weights", required=FALSE),
        penalty = penalty,
        scale = scale
    )

    processed <- hardhat::mold(form_combined, data, blueprint=bp)
    ei_riesz_bridge(processed, ...)
}



#' @export
#' @rdname ei_riesz
ei_riesz.ei_spec <- function(x, weights, penalty, scale=TRUE, ...) {
    spec = x
    validate_ei_spec(spec)

    form = as.formula(paste0(
        paste0(attr(spec, "ei_y"), collapse=" + "), " ~ ",
        paste0(c(attr(spec, "ei_x"), attr(spec, "ei_z")), collapse=" + ")
    ))

    bp = hardhat::new_default_formula_blueprint(
        intercept = FALSE,
        composition = "matrix",
        indicators = "one_hot",
        ei_x = attr(spec, "ei_x"),
        ei_n = attr(spec, "ei_n"),
        ei_wgt = check_make_weights(!!enquo(weights), spec, arg="weights", required=FALSE),
        penalty = penalty,
        scale = scale
    )

    processed <- hardhat::mold(form, spec, blueprint=bp)
    ei_riesz_bridge(processed, ...)
}


#' @export
#' @rdname ei_riesz
ei_riesz.data.frame <- function(x, z, total, weights, penalty, scale=TRUE, ...) {
    if (length(both <- intersect(colnames(x), colnames(z))) > 0) {
        cli_abort(c("Predictors and covariates must be distinct",
                    ">"="Got: {.var {both}}"), call=parent.frame())
    }
    if (!(is.matrix(z) || is.data.frame(z))) {
        cli_abort("{.arg z} must be a matrix or data frame.", call=parent.frame())
    }

    bp = hardhat::new_default_xy_blueprint(
        intercept = FALSE,
        composition = "matrix",
        ei_x = colnames(x),
        ei_n = check_make_weights(!!enquo(total)),
        ei_wgt = check_make_weights(!!enquo(weights), arg="weights", required=FALSE),
        penalty = penalty,
        scale = scale
    )
    x = cbind(x, z)

    processed <- hardhat::mold(x, NULL, blueprint=bp)
    ei_riesz_bridge(processed, ...)
}

#' @export
#' @rdname ei_riesz
ei_riesz.matrix <- function(x, z, total, weights, penalty, scale=TRUE, ...) {
    ei_riesz.data.frame(x, z, total, weights, penalty, scale, ...)
}


#' @export
#' @rdname ei_riesz
ei_riesz.default <- function(x, ...) {
    if (missing(x))
        cli_abort("{.fn ei_riesz} requires arguments.", call=NULL)
    cli_abort("{.fn ei_riesz} is not defined for a {.cls {class(x)}}.", call=NULL)
}


# Bridge and implementation ---------------------------------------------------

ei_riesz_bridge <- function(processed, ...) {
    err_call = rlang::new_call(rlang::sym("ei_riesz"))
    xz = processed$predictors
    idx_x = match(processed$blueprint$ei_x, colnames(xz))
    z = xz[, -idx_x, drop=FALSE]
    x = pull_x(xz, idx_x)
    check_preds(x, call=err_call)
    total = processed$blueprint$ei_n
    weights = processed$blueprint$ei_wgt
    penalty = processed$blueprint$penalty

    # normalize
    z_shift = colSums(z * weights) / sum(weights)
    z = shift_cols(z, z_shift)
    if (isTRUE(processed$blueprint$scale)) {
        z_scale = (colSums(z^2 * weights) / sum(weights))^-0.5
        z = scale_cols(z, z_scale)
    } else {
        z_scale = rep(1, ncol(z))
    }

    # NA checking
    if (any(is.na(x))) cli_abort("Missing values found in predictors.", call=err_call)
    if (any(is.na(z))) cli_abort("Missing values found in covariates.", call=err_call)

    fit <- ei_riesz_impl(x, z, total, weights, penalty)

    new_ei_riesz(
        weights = fit$alpha,
        weights_loo = fit$loo,
        nu2 = fit$nu2,
        penalty = penalty,
        z_shift = z_shift,
        z_scale = z_scale,
        blueprint = processed$blueprint
    )
}

#' @param total A vector of total observations per unit.
#'
#' @rdname ei-impl
#' @export
ei_riesz_impl <- function(x, z, total, weights=rep(1, nrow(x)), penalty) {
    int_scale = 1 + 1e2*sqrt(penalty)
    w = weights / mean(weights)
    xz = row_kronecker(x, z, int_scale)
    sqrt_w = sqrt(weights / mean(weights))
    udv = svd(xz * sqrt_w)

    alpha = matrix(nrow=nrow(x), ncol=ncol(x))
    loo = matrix(nrow=nrow(x), ncol=ncol(x))
    nu2 = numeric(ncol(x))
    for (group in seq_len(ncol(x))) {
        fit = riesz_svd(xz, udv, ncol(z), total, w, sqrt_w, group, penalty)
        alpha[, group] = fit$alpha * int_scale * w
        loo[, group] = fit$loo * int_scale * w
        nu2[group] = fit$nu2  * int_scale^2
    }
    colnames(alpha) = colnames(x)
    colnames(loo) = colnames(x)
    names(nu2) = colnames(x)

    list(alpha = alpha, loo = loo, nu2 = nu2)
}

# Model type ------------------------------------------------------------------

new_ei_riesz <- function(..., blueprint) {
    hardhat::new_model(..., blueprint = blueprint, class = "ei_riesz")
}

#' @export
print.ei_riesz <- function(x, ...) {
    cat_line("A Riesz representer with ",
             ncol(x$weights), " groups and ",
             nrow(x$weights), " observations")
    cat_line("Fit with penalty = ", signif(x$penalty))
}

#' @export
summary.ei_riesz <- function(object, ...) {
    cat_line("Second moment of representer:")
    print(colMeans(object$weights^2))
    cat_line()
    cat_line("Second moment of representer (leave-one-out):")
    print(colMeans(object$weights_loo^2))
}

#' Extract Riesz representer weights
#'
#' Extracts a single set of Riesz representer weights from an `ei_riesz` object,
#' for a selected group.
#'
#' @param object An [ei_riesz()] object.
#' @param group The group for which to extract the weights, as a numeric index
#'   or a character column name. The special (default) value `TRUE` will return
#'   a matrix of weights, with each column corresponding to a group.
#' @param loo If `TRUE`, return the leave-one-out weights
#' @param ... Additional arguments (ignored)
#'
#' @returns A numeric vector of weights
#' @export
weights.ei_riesz <- function(object, group=TRUE, loo=FALSE, ...) {
    if (isTRUE(group)) {
        group = seq_len(ncol(object$weights))
    }

    if (isTRUE(loo)) {
        object$weights_loo[, group]
    } else {
        object$weights[, group]
    }
}

