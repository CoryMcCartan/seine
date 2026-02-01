# User-facing model functions --------------------------------------------------

#' Fit an ecological inference regression model
#'
#' Fits a penalized regression model for ecological inference, allowing for
#' overall and unit-level estimates of conditional means using [ei_est()].
#'
#' The regression is calculated using the singular value decomposition, which
#' allows for efficient recalculation under different `penalty` values as part
#' of leave-one-out cross-validation.
#' When `bounds` are provided, the regression is calculated via quadratic
#' programming, as there is no closed-form solution. The unbounded regression
#' is run to select the `penalty` automatically in this case, if it is not
#' provided. Estimation is still efficient, though somewhat slower than in the
#' unbounded case. The covariance matrix of the estimates is not available when
#' bounds are applied.
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
#'   Covariates in an [ei_spec] object are shifted to have mean zero. If
#'   `scale=TRUE` (the default), they are also scaled to have unit variance.
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
#'   These are shifted to have mean zero. If `scale=TRUE` (the default), they
#'   are also scaled to have unit variance.
#' @param data When a **formula** is used, `data` is a **data frame** containing
#'   both the predictors and the outcome.
#' @param formula A formula such as `y ~ x0 + x1 | z` specifying the outcome `y`
#'  regressed on the predictors of interest `x` and any covariates `z`.
#'  The predictors should form a partition, that is, `x0 + x1 = 1` for each
#'  observation. Users can be include more than two predictors as well, e.g.
#'  `pct_white + pct_black + pct_hisp + pct_other`.
#'  If there are just two predictors, it is acceptable to only include one in
#'  the formula; the other will be formed as 1 minus the provided predictor.
#'  Include additional covariates separated by a vertical bar `|`.
#'  These covariates are strongly recommended for reliable ecological inference.
#'  Covariates are shifted to have mean zero. If `scale=TRUE` (the default),
#'  they are also scaled to have unit variance.
#' @param weights <[`data-masking`][rlang::args_data_masking]> A vector of unit
#'   weights for estimation. These may be the same or different from the total
#'   number of observations in each aggregate unit (see the `total` argument to
#'   [ei_spec()]). See the discussion below under 'Weights' for choosing this
#'   parameter. The default, uniform weights, makes a slightly
#'   stronger-than-necessary assumption about the relationship between the
#'   unit totals and the unknown data.
#' @param penalty The ridge penalty (a non-negative scalar). Set to `NULL` to
#'    automatically determine the penalty which minimizes mean-square error,
#'    via an efficient leave-one-out cross validation procedure.
#'    The ridge regression solution is
#'    \deqn{\hat\beta = (X^\top X + \lambda I)^{-1}X^\top y,}
#'    where \eqn{\lambda} is the value of `penalty`.
#'    One can equivalently think of the penalty as imposing a
#'    \eqn{\mathcal{N}(0, \sigma^2/\lambda^2)} prior on the \eqn{\beta}.
#'    Keep in mind when choosing `penalty` manually that covariates in `z` are
#'    scaled to have mean zero and unit variance before fitting.
#' @param bounds A vector `c(min, max)` of bounds for the outcome.
#'   If `bounds = NULL`, they will be inferred from the outcome variable:
#'   if it is contained within \eqn{[0, 1]}, for instance, then the bounds will
#'   be `c(0, 1)`. The default `bounds = FALSE` uses an unbounded outcome.
#' @param sum_one If `TRUE`, the outcome variables are constrained to sum to one.
#'   Can only apply when `bounds` are enforced and there is more than one
#'   outcome variable. If `NULL`, infers `sum_one = TRUE` when the bounds
#'   are `c(0, 1)` the outcome variables sum to 1.
#' @param scale If `TRUE`, scale covariates `z` to have unit variance.
#' @param vcov If `TRUE`, calculate and return the a scaled covariance matrix of
#'    the estimated coefficients. When `bounds` are provided, the (scaled)
#'    covariance matrix for the unbounded estimate is returned as a conservative
#'    approximation.
#'    The covariance matrix is "scaled" because it does not include the
#'    residual variance. For the covariance for a particular outcome variable,
#'    multiply the returned `$vcov_u` by `sigma2` for that outcome.
#' @param ... Not currently used, but required for extensibility.
#'
#' @returns An `ei_ridge` object, which supports various [ridge-methods].
#'
#' @inherit ei_est references
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
#'                total = pres_total, covariates = c(pop_urban, farm))
#' ei_ridge(spec)
#'
#' ei_ridge(pres_dem_hum + pres_rep_nix + pres_ind_wal + pres_abs ~
#'       vap_white + vap_black + vap_other | pop_urban + farm, data = elec_1968)
#'
#' # bounds inferred
#' all.equal(
#'   fitted(ei_ridge(spec, bounds = NULL)),
#'   fitted(ei_ridge(spec, bounds = 0:1))
#' )
#'
#' # bounds enforced
#' min(fitted(ei_ridge(spec)))
#' min(fitted(ei_ridge(spec, bounds = 0:1)))
#' @export
ei_ridge <- function(x, ..., weights, bounds = FALSE, sum_one = FALSE, penalty = NULL, scale = TRUE, vcov = TRUE) {
    UseMethod("ei_ridge")
}


#' @export
#' @rdname ei_ridge
ei_ridge.formula <- function(formula, data, weights, bounds=FALSE, sum_one=FALSE,
                             penalty=NULL, scale=TRUE, vcov=TRUE, ...) {
    forms = ei_forms(formula)
    form_preds = terms(rlang::new_formula(lhs=NULL, rhs=forms$predictors))
    form_combined = rlang::new_formula(forms$outcome, expr(!!forms$predictors + !!forms$covariates))

    bp = hardhat::new_default_formula_blueprint(
        intercept = FALSE,
        composition = "matrix",
        indicators = "one_hot",
        ei_x = attr(form_preds, "term.labels"),
        ei_wgt = check_make_weights(!!enquo(weights), data, arg="weights", required=FALSE),
        bounds = bounds,
        sum_one = sum_one,
        penalty = penalty,
        scale = scale,
        subclass = "ei_ridge_blueprint"
    )

    processed <- hardhat::mold(form_combined, data, blueprint=bp)
    ei_ridge_bridge(processed, vcov, ...)
}



#' @export
#' @rdname ei_ridge
ei_ridge.ei_spec <- function(x, weights, bounds=FALSE, sum_one=FALSE, penalty=NULL,
                             scale=TRUE, vcov=TRUE, ...) {
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
        ei_wgt = check_make_weights(!!enquo(weights), spec, arg="weights", required=FALSE),
        bounds = bounds,
        sum_one = sum_one,
        penalty = penalty,
        scale = scale,
        subclass = "ei_ridge_blueprint"
    )

    processed <- hardhat::mold(form, spec, blueprint=bp)
    ei_ridge_bridge(processed, vcov, ...)
}


#' @export
#' @rdname ei_ridge
ei_ridge.data.frame <- function(x, y, z, weights, bounds=FALSE, sum_one=FALSE, penalty=NULL,
                                scale=TRUE, vcov=TRUE, ...) {
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
        bounds = bounds,
        sum_one = sum_one,
        penalty = penalty,
        scale = scale,
        subclass = "ei_ridge_blueprint"
    )
    x = cbind(x, z)

    processed <- hardhat::mold(x, y, blueprint=bp)
    ei_ridge_bridge(processed, vcov, ...)
}

#' @export
#' @rdname ei_ridge
ei_ridge.matrix <- function(x, y, z, weights, bounds=FALSE, sum_one=FALSE, penalty=NULL,
                            scale=TRUE, vcov=TRUE, ...) {
    ei_ridge.data.frame(x, y, z, weights, penalty, sum_one, bounds, scale, vcov, ...)
}


#' @export
#' @rdname ei_ridge
ei_ridge.default <- function(x, ...) {
    if (missing(x))
        cli_abort("{.fn ei_ridge} requires arguments.", call=NULL)
    cli_abort("{.fn ei_ridge} is not defined for a {.cls {class(x)}}.", call=NULL)
}

# Bridge and implementation ---------------------------------------------------


# Creates bounds _after_ molding outcomes
#' @exportS3Method hardhat::run_mold
run_mold.ei_ridge_blueprint <- function(blueprint, ...) {
    processed = NextMethod("run_mold")

    # update bounds
    bounds = ei_bounds(processed$blueprint$bounds, processed$outcomes)
    processed$blueprint = hardhat::update_blueprint(
        processed$blueprint,
        bounds = bounds
    )

    processed
}

ei_ridge_bridge <- function(processed, vcov, ...) {
    err_call = rlang::new_call(rlang::sym("ei_ridge"))
    bp = processed$blueprint
    xz = processed$predictors
    idx_x = match(bp$ei_x, colnames(xz))
    x = pull_x(xz, idx_x)
    check_preds(x, call=err_call)
    z = xz[, -idx_x, drop=FALSE]
    weights = bp$ei_wgt

    # normalize
    z_shift = colSums(z * weights) / sum(weights)
    z = shift_cols(z, z_shift)
    if (isTRUE(bp$scale)) {
        z_scale = (colSums(z^2 * weights) / sum(weights))^-0.5
        z = scale_cols(z, z_scale)
    } else {
        z_scale = rep(1, ncol(z))
    }

    y = as.matrix(processed$outcomes)

    # NA checking
    if (any(is.na(x))) cli_abort("Missing values found in predictors.", call=err_call)
    if (any(is.na(y))) cli_abort("Missing values found in outcome.", call=err_call)
    if (any(is.na(z))) cli_abort("Missing values found in covariates.", call=err_call)

    if (ncol(z) == 0) {
        bp$penalty = 0
    }
    if (is.null(bp$sum_one) && all(bp$bounds == c(0, 1))) {
        bp$sum_one = isTRUE(all.equal(rowSums(y), rep(1, nrow(y))))
    }

    fit <- ei_ridge_impl(x, y, z, weights, bp$bounds, bp$sum_one, bp$penalty, vcov)

    new_ei_ridge(
      coef = fit$coef,
      y = y,
      fitted = fit$fitted,
      vcov_u = fit$vcov_u,
      sigma2 = fit$sigma2,
      r2 = diag(as.matrix(cor(fit$fitted, y)^2)),
      penalty = fit$penalty,
      int_scale = fit$int_scale,
      z_shift = z_shift,
      z_scale = z_scale,
      blueprint = bp
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
#' @param bounds A vector `c(min, max)` of bounds for the outcome.
#' @param penalty The ridge penalty (a non-negative scalar), which must be
#'   specified for [ei_riesz_impl()] but can be automatically estimated with
#'   [ei_ridge_impl()] by providing `penalty=NULL`.
#' @inheritParams ei_ridge
#'
#' @returns A list with model components.
#'
#' @rdname ei-impl
#' @export
ei_ridge_impl <- function(x, y, z, weights=rep(1, nrow(x)),
                          bounds=c(-Inf, Inf), sum_one=NULL, penalty=NULL, vcov=TRUE) {
    int_scale = if (!is.null(penalty) && penalty == 0) 1 + 1e2*sqrt(penalty) else 1e4
    xz = row_kronecker(x, z, int_scale)
    sqrt_w = sqrt(weights / mean(weights))
    udv = svd(xz * sqrt_w)

    vcov = isTRUE(vcov)
    enforce = is.finite(bounds)
    if (!any(enforce)) { # unbounded
        if (isTRUE(sum_one)) {
            cli_abort("{.fn ei_ridge} cannot enforce sum-to-one constraint when outcome is unbounded.")
        }
        fit = if (is.null(penalty)) {
            ridge_auto(udv, y, sqrt_w, vcov)
        } else {
            ridge_svd(udv, y, sqrt_w, penalty, vcov)
        }
    } else {
        if (is.null(penalty) || vcov) {
            unb_fit = ridge_auto(udv, y, sqrt_w, vcov)
        }
        if (is.null(penalty)) {
            penalty = unb_fit$penalty
        }

        fit = ridge_bounds(xz, z, y, weights, bounds, sum_one, penalty)
        if (vcov) {
            fit$vcov_u = unb_fit$vcov_u
        }
    }

    rownames(fit$coef) = colnames(xz)
    if (!is.null(fit$vcov_u)) {
        rownames(fit$vcov_u) = colnames(fit$vcov_u) = colnames(xz)
    }
    names(fit$sigma2) = colnames(y)
    fit$int_scale = int_scale

    fit
}

predict_ei_ridge_bridge <- function(type, object, processed) {
    type = rlang::arg_match(type, "numeric")

    x = processed$predictors
    idx_x = match(object$blueprint$ei_x, colnames(x))
    z = x[, -idx_x, drop=FALSE]
    x = pull_x(x, idx_x)
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

# helper to pull columns from predictor matrix and add .other if needed
# used by ei_riesz and ei_est as well
pull_x <- function(x, idx_x) {
    if (any(is.na(idx_x))) {
        cli_abort(c("Specification error.",
                    "i"="Factors are not allowed as predictors.
                          Check your formula or specification.",
                    ">"="If there are no errors in your specification, please report
                         this at {.url https://github.com/CoryMcCartan/seine/issues}"),
                  call = NULL)
    }

    x = x[, idx_x, drop=FALSE]
    if (ncol(x) == 1) {
        x = cbind(x, 1 - x)
        colnames(x)[2] = ".other"
    }
    x
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
    bounds = x$blueprint$bounds
    if (any(is.finite(bounds))) {
        sumt1 = if (isTRUE(x$blueprint$sum_one)) " and constrained to sum to 1" else ""
        pl = if (ncol(x$y) > 1) "s" else ""
        cat_line("With outcome", pl, " bounded in (", bounds[1], ", ", bounds[2], ")", sumt1)
    }
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

#' @describeIn ridge-methods Extract covariance of coefficient estimates.
#' @export
vcov.ei_ridge <- function(object, ...) {
    object$vcov_u
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
