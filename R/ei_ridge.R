# User-facing model functions --------------------------------------------------

#' Fit an ecological inference regression model
#'
#' Fits a penalized regression model for ecological inference, allowing for
#' overall and unit-level estimates of conditional means using [ei_est()].
#'
#' The regression is calculated using the singular value decomposition, which
#' allows for efficient recalculation under different `penalty` values as part
#' of leave-one-out cross-validation.
#'
#' @section Weights:
#' The weakest identification result for ecological inference makes no
#' assumption about the number of observations per aggregate unit (the totals).
#' It requires, however, weighting the estimation steps according to the totals.
#' This may reduce efficiency when the totals are variable and a slightly
#' stronger condition holds.
#'
#' Specifically, if the totals are conditionally mean-independent of the missing
#' data (the aggregation-unit level means of the outcome within each predictor
#' level), given covariates, then it is appropriate to use uniform weights in
#' estimation, or any fixed set of weights.
#'
#' In general, estimation efficiency is improved when units with larger variance
#' in the outcome receive less weight. Various bulit-in options are provided by
#' the helper functions in [ei_wgt()].
#'
#'
#' @param x Depending on the context:
#'   * A **data frame** of predictors.
#'   * A **matrix** of predictors.
#'   * An [ei_spec] object containing the outcome, predictor, and covariates.
#'
#'   Predictors must be proportions that sum to 1 across rows.
#'   You can use [ei_proportions()] to assist in preparing predictor variables.
#' @param y When `x` is a **data frame** or **matrix**, `y` is the outcome
#' specified as:
#'   * A **data frame** with numeric columns.
#'   * A **matrix**
#'   * A numeric **vector**.
#'
#'   When the outcome is a proportion, you can use [ei_proportions()] to assist
#'   in preparing it.
#' @param z When `x` is a **data frame** or **matrix**, `w` are any covariates,
#' specified as:
#'   * A **data frame** with numeric columns.
#'   * A **matrix**
#'
#'   These are scaled internally to have mean zero and unit variance.
#' @param data When a **formula** is used, `data` is a **data frame** containing
#'   both the predictors and the outcome.
#' @param formula A formula such as `y | x ~ z` specifying the outcome and
#'   predictor terms on the left-hand side, and any covariate terms on the
#'   right-hand side.
#'   The outcome and predictor variables must be separated by a vertical bar
#'   `|` on the left-hand side.
#'   This is because EI models the conditional mean of the outcome given
#'   predictors, and how this may vary with covariates.
#'   See the examples for more details.
#' @param weights <[`data-masking`][rlang::args_data_masking]> A vector of unit
#'   weights for estimation. These may be the same or different from the total
#'   number of observations in each aggregate unit (see the `total` argument to
#'   [ei_spec()]). See the discussion below under 'Weights' for choosing this
#'   parameter. The default, uniform weights, makes a slightly
#'   stronger-than-necessary assumption about the relationship between the
#'   unit totals and the unknown data.
#' @param penalty The ridge penalty. Set to `NULL` to automatically determine
#'    the penalty which minimizes mean-square error, via an efficient
#'    leave-one-out cross validation procedure. The ridge regression solution is
#'    \deqn{\hat\beta = (X^\top X + \lambda I)^{-1}X^\top y,}
#'    where \eqn{\lambda} is the value of `penalty`.
#'    Keep in mind when choosing `penalty` manually that covariates in `z` are
#'    scaled to have mean zero and unit variance before fitting.
#' @param ... Not currently used, but required for extensibility.
#'
#' @returns An `ei_ridge` object, which supports various [ridge-methods].
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
#'                total = pres_total, covariates = c(pop_urban, farm))
#' ei_ridge(spec)
#'
#' ei_ridge(pres_dem_hum + pres_rep_nix + pres_ind_wal + pres_abs + pres_oth ~
#'       vap_white + vap_black + vap_other | pop_urban + farm, data = elec_1968)
#'
#' @export
ei_ridge <- function(x, ...) {
    UseMethod("ei_ridge")
}



#' @export
#' @rdname ei_ridge
ei_ridge.formula <- function(formula, data, weights, penalty=NULL, ...) {
    forms = ei_forms(formula)
    form_preds = terms(rlang::new_formula(lhs=NULL, rhs=forms$predictors))
    form_combined = rlang::new_formula(forms$outcome, expr(!!forms$predictors + !!forms$covariates))

    bp = hardhat::new_default_formula_blueprint(
        intercept = FALSE,
        composition = "matrix",
        ei_x = attr(form_preds, "term.labels"),
        ei_wgt = check_make_weights(!!enquo(weights), data, arg="weights", required=FALSE),
        penalty = penalty,
    )

    processed <- hardhat::mold(form_combined, data, blueprint=bp)
    ei_ridge_bridge(processed, ...)
}



#' @export
#' @rdname ei_ridge
ei_ridge.ei_spec <- function(x, weights, penalty=NULL, ...) {
    spec = x
    x = spec[c(attr(spec, "ei_x"), attr(spec, "ei_z"))]
    # handle factors
    chr_cols = logical(length(x))
    for (col in seq_len(length(x))) {
        if (is.character(x[[col]]) || is.factor(x[[col]])) {
            chr_cols[col] = TRUE
            form = formula(paste("~ 0 +", names(x)[col]))
            x = cbind(x, model.matrix(form, x))
        }
    }
    x = x[, !chr_cols, drop=FALSE]

    y = spec[attr(spec, "ei_y")]

    bp = hardhat::new_default_xy_blueprint(
        intercept = FALSE,
        composition = "matrix",
        ei_x = attr(spec, "ei_x"),
        ei_wgt = check_make_weights(!!enquo(weights), x, arg="weights", required=FALSE),
        penalty = penalty,
    )

    processed = hardhat::mold(x, y, blueprint=bp)
    ei_ridge_bridge(processed, ...)
}


#' @export
#' @rdname ei_ridge
ei_ridge.data.frame <- function(x, y, z, weights, penalty=NULL, ...) {
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
        ei_wgt = check_make_weights(!!enquo(weights), arg="weights", required=FALSE),
        penalty = penalty,
    )
    x = cbind(x, z)

    processed <- hardhat::mold(x, y, blueprint=bp)
    ei_ridge_bridge(processed, ...)
}

#' @export
#' @rdname ei_ridge
ei_ridge.matrix <- function(x, y, z, weights, penalty=NULL, ...) {
    ei_ridge.data.frame(x, y, z, weights, penalty, ...)
}


#' @export
#' @rdname ei_ridge
ei_ridge.default <- function(x, ...) {
    if (missing(x))
        cli_abort("{.fn ei_ridge} requires arguments.", call=NULL)
    cli_abort("{.fn ei_ridge} is not defined for a {.cls {class(x)}}.", call=NULL)
}

# Bridge and implementation ---------------------------------------------------

ei_ridge_bridge <- function(processed, ...) {
    x = processed$predictors
    idx_x = match(processed$blueprint$ei_x, colnames(x))
    z = x[, -idx_x, drop=FALSE]
    x = x[, idx_x, drop=FALSE]
    if (ncol(x) == 1) {
        x = cbind(x, 1 - x)
        colnames(x)[2] = ".other"
    }
    weights = processed$blueprint$ei_wgt

    # normalize
    z_shift = colSums(z * weights) / sum(weights)
    z = shift_cols(z, z_shift)
    z_scale = (colSums(z^2 * weights) / sum(weights))^-0.5
    z = scale_cols(z, z_scale)

    y = as.matrix(processed$outcomes)

    # NA checking
    if (any(is.na(x))) cli_abort("Missing values found in predictors.")
    if (any(is.na(y))) cli_abort("Missing values found in outcome.")
    if (any(is.na(z))) cli_abort("Missing values found in covariates.")

    if (ncol(z) == 0)
        processed$blueprint$penalty = 0

    fit <- ei_ridge_impl(x, y, z, weights, processed$blueprint$penalty)

    new_ei_ridge(
      coef = fit$coef,
      y = y,
      fitted = fit$fitted,
      r2 = diag(as.matrix(cor(fit$fitted, y)^2)),
      penalty = fit$penalty,
      int_scale = fit$int_scale,
      z_shift = z_shift,
      z_scale = z_scale,
      blueprint = processed$blueprint
    )
}

#' Low-level implementations of `ei_ridge()` and `ei_riesz()`
#'
#' No checks are performed on the inputs.
#' Use of [ei_ridge()] and [ei_riesz()] is strongly recommended unless many
#' regressions must be fit, e.g., within a tight loop.
#' Only works for a single outcome, i.e., `y` must be a vector, not a matrix.
#'
#' @param x A matrix of predictors
#' @param y A vector of outcomes
#' @param z A matrix of covariates
#' @param weights A vector of estimation weights
#' @inheritParams ei_riesz
#'
#' @returns A list with model components.
#'
#' @rdname ei-impl
#' @export
ei_ridge_impl <- function(x, y, z, weights, penalty=NULL) {
    int_scale = if (!is.null(penalty) && penalty == 0) 1 + 1e2*sqrt(penalty) else 1e4
    xz = row_kronecker(x, z, int_scale)
    sqrt_w = sqrt(weights / mean(weights))
    udv = svd(xz * sqrt_w)

    fit = if (is.null(penalty)) {
        ridge_auto(udv, y, sqrt_w)
    } else {
        ridge_svd(udv, y, sqrt_w, penalty)
    }

    rownames(fit$coef) = colnames(xz)
    fit$int_scale = int_scale

    fit
}

predict_ei_ridge_bridge <- function(type, object, processed) {
    type = rlang::arg_match(type, "numeric")

    x = processed$predictors
    idx_x = match(object$blueprint$ei_x, colnames(x))
    z = x[, -idx_x, drop=FALSE]
    x = x[, idx_x, drop=FALSE]
    if (ncol(x) == 1) {
        x = cbind(x, 1 - x)
        colnames(x)[2] = ".other"
    }

    # normalize
    z = shift_cols(z, object$z_shift)
    z = scale_cols(z, object$z_scale)

    # NA checking
    if (any(is.na(x))) cli_abort("Missing values found in predictors.")
    if (any(is.na(z))) cli_abort("Missing values found in covariates.")

    switch(
        type,
        numeric = predict_ei_ridge_numeric(object, x, z)
    )
}

predict_ei_ridge_numeric <- function(object, x, z) {
    xz = row_kronecker(x, z, object$int_scale)
    pred = as.list(as.data.frame(xz %*% object$coef))
    do.call(hardhat::spruce_numeric_multiple, pred)
}

# Model type ------------------------------------------------------------------

new_ei_ridge <- function(..., blueprint) {
    hardhat::new_model(..., blueprint = blueprint, class = "ei_ridge")
}

#' @export
print.ei_ridge <- function(x, ...) {
    cat_line("An ecological inference model with ",
             ncol(x$coef), " outcomes, ",
             length(x$blueprint$ei_x), " groups, and ",
             nrow(x$fitted), " observations")
    cat_line("Fit with penalty = ", signif(x$penalty))
}

#' Methods for [`ei_ridge`] models
#'
#' Models fitted with [ei_ridge()] support various generic methods.
#'
#' @param object A fitted [ei_ridge] model
#' @param ... Additional arguments (ignored)
#'
#' @name ridge-methods
NULL

#' @describeIn ridge-methods Predict from an `ei_ridge` model.
#' @param new_data A data frame, matrix, or [ei_spec] of new predictors.
#' @param type The type of predictions to generate; only `"numeric"` is supported.
#' @export
predict.ei_ridge <- function(object, new_data, type="numeric", ...) {
    processed = hardhat::forge(new_data, object$blueprint)
    out = predict_ei_ridge_bridge(type, object, processed)
    hardhat::validate_prediction_size(out, new_data)
    out
}


#' @describeIn ridge-methods Extract fitted values.
#' @export
fitted.ei_ridge <- function(object, ...) {
    object$fitted
}

#' @describeIn ridge-methods Extract residuals.
#' @export
residuals.ei_ridge <- function(object, ...) {
    object$y - object$fitted
}


#' @describeIn ridge-methods Summarize the model's fitted values and \eqn{R^2}.
#' @export
summary.ei_ridge <- function(object, ...) {
    cat_line("Fitted values:")
    print(summary((object$fitted)))
    cat_line()
    cat_line("R-squared by outcome:")
    print(object$r2)
}

#' @describeIn ridge-methods Extract estimation weights from a fitted model.
#' @param normalize If `TRUE`, normalize the weights to have mean 1.
#' @export
weights.ei_ridge <- function(object, normalize = TRUE, ...) {
    w = object$blueprint$ei_wgt
    if (is.null(w)) {
        w = rep(1, nrow(object))
    }
    if (isTRUE(normalize)) {
        w / mean(w)
    } else{
        w
    }
}
