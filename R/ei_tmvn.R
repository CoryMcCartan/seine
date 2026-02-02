#' Fit King's ecological inference model
#'
#' Fits King's (1997) regression model for ecological inference, allowing for overall and
#' unit-level estimates of conditional means using [ei_est()].
#'
#' @inheritParams ei_ridge
#' @param y When `x` is a **data frame** or **matrix**, `y` is the outcome
#' specified as:
#'   * A **data frame** with numeric columns.
#'   * A **matrix**
#'   * A numeric **vector**.
#'
#'   When the outcome is a proportion, you can use [ei_proportions()] to assist
#'   in preparing it.
#'   For functions which incorporate bounds on the outcome variable, these are
#'   specified by the `bounds` parameter, which is inferred from context, by
#'   default.
#' @param bounds A vector `c(min, max)` of bounds for the outcome. The default
#'   `NULL` will use the strictest bounds consistent with the observed data:
#'   between 0 and 1, where possible, then greater than 0, then unbounded.
#' @param ... Not currently used, but required for extensibility.
#'
#' @returns An `ei_tmvn` object.
#'
#' @references
#' King, G. (1997). *A solution to the ecological inference problem:
#' Reconstructing individual behavior from aggregate data.* Princeton University
#' Press.
#'
#' @examples
#' data(elec_1968)
#'
#' \dontrun{
#' ei_tmvn(pres_ind_wal ~ vap_black, elec_1968, weights=pres_total)
#' }
#' @noRd
ei_tmvn <- function(x, ...) {
    UseMethod("ei_tmvn")
}


# @export
# @rdname ei_tmvn
#' @export
ei_tmvn.formula <- function(formula, data, weights, bounds=NULL, penalty=0, scale=TRUE, ...) {
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
        penalty = penalty,
        scale = scale,
        subclass = "ei_tmvn_blueprint"
    )

    processed <- hardhat::mold(form_combined, data, blueprint=bp)
    ei_tmvn_bridge(processed, ...)
}


# @export
# @rdname ei_tmvn
#' @export
ei_tmvn.ei_spec <- function(x, weights, bounds=NULL, penalty=0, scale=TRUE, ...) {
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
        penalty = penalty,
        scale = scale,
        subclass = "ei_tmvn_blueprint"
    )

    processed <- hardhat::mold(form, spec, blueprint=bp)
    ei_tmvn_bridge(processed, ...)
}


# @export
# @rdname ei_tmvn
#' @export
ei_tmvn.data.frame <- function(x, y, z, weights, bounds=NULL, penalty=0, scale=TRUE, ...) {
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
        penalty = penalty,
        scale = scale,
        subclass = "ei_tmvn_blueprint"
    )
    x = cbind(x, z)

    processed <- hardhat::mold(x, y, blueprint=bp)
    ei_tmvn_bridge(processed, ...)
}

# @export
# @rdname ei_tmvn
#' @export
ei_tmvn.matrix <- function(x, y, z, weights, bounds=NULL, penalty=0, scale=TRUE, ...) {
    ei_tmvn.data.frame(x, y, z, weights, bounds, penalty, scale, ...)
}


# @export
# @rdname ei_tmvn
#' @export
ei_tmvn.default <- function(x, ...) {
    if (missing(x))
        cli_abort("{.fn ei_tmvn} requires arguments.", call=NULL)
    cli_abort("{.fn ei_tmvn} is not defined for a {.cls {class(x)}}.", call=NULL)
}

# Bridge ----------------------------------------------------------------------

# Creates bounds _after_ molding outcomes
#' @exportS3Method hardhat::run_mold
run_mold.ei_tmvn_blueprint <- function(blueprint, ...) {
    processed = NextMethod("run_mold")

    # update bounds
    bounds = check_bounds(processed$blueprint$bounds, processed$outcomes)
    if (bounds[1] == -Inf && bounds[2] == Inf) {
        cli_abort(c("Bounds were set or inferred to be `c(-Inf, Inf)`.",
                    "i"="Use {.fn ei_ridge} for unbounded regression."),
                  call=rlang::new_call(rlang::sym("ei_tmvn")))
    }

    processed$blueprint = hardhat::update_blueprint(
        processed$blueprint,
        bounds = bounds
    )

    processed
}

ei_tmvn_bridge <- function(processed, ...) {
    err_call = rlang::new_call(rlang::sym("ei_tmvn"))
    bp = processed$blueprint
    xz = processed$predictors
    idx_x = match(bp$ei_x, colnames(xz))
    z = xz[, -idx_x, drop=FALSE]
    x = pull_x(xz, idx_x)
    check_preds(x, call=err_call)
    weights = processed$blueprint$ei_wgt

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
    # standardize bounds
    if (all(is.finite(bp$bounds))) {
        y_shift = bp$bounds[1]
        y_scale = bp$bounds[2] - bp$bounds[1]
        y = (y - y_shift) / y_scale
    } else {
        cli_abort(c("Inference on half-open regions is not currently supported.",
                    "x"="Bounds must be finite; got {.val {bp$bounds}}.",
                    "i"="{.fn ei_ridge} supports half-open bounds."),
                     call=err_call)
    }

    # NA checking
    if (any(is.na(x))) cli_abort("Missing values found in predictors.", call=err_call)
    if (any(is.na(y))) cli_abort("Missing values found in outcome.", call=err_call)
    if (any(is.na(z))) cli_abort("Missing values found in covariates.", call=err_call)

    if (ncol(z) == 0) {
        bp$penalty = 0
    }

    fit <- ei_tmvn_impl(x, y, z, weights, bp$bounds, bp$penalty)

    do.call(new_ei_tmvn, c(fit, list(
        y = y,
        penalty = bp$penalty,
        y_shift = y_shift,
        y_scale = y_scale,
        z_shift = z_shift,
        z_scale = z_scale,
        blueprint = bp
    )))
    # new_ei_tmvn(
    #   coefs = fit$coefs,
    #   blueprint = processed$blueprint
    # )
}

# Model type ------------------------------------------------------------------

new_ei_tmvn <- function(..., blueprint) {
    hardhat::new_model(..., blueprint = blueprint, class = "ei_tmvn")
}

#' @export
print.ei_tmvn <- function(x, ...) {
    cat_line("An ecological inference model")
    cat_line("<placeholder print method>")
    if ("b_global" %in% names(x)) {
        cat_line("Global estimate:")
        print(x$b_global)
    } else {
        cat_line("Global estimate (no draws):")
        print(x$est$beta[, 1])
    }
}
