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
#' @param bounds The bounds to enforce on the outcome, as a length-two vector
#'   `c(lower, upper)`. The default `NULL` will use the strictest bounds
#'   consistent with the observed data: between 0 and 1, where possible, then
#'   greater than 0, then unbounded.
#' @param ... Not currently used, but required for extensibility.
#'
#' @returns An `ei_model` object.
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
#' ei_model(pres_ind_wal ~ vap_white + vap_black, elec_1968, weights=pres_total)
#' }
#' @export
ei_model <- function(x, ...) {
  UseMethod("ei_model")
}


#' @export
#' @rdname ei_model
ei_model.formula <- function(formula, data, weights, bounds=NULL, ...) {
    forms = ei_forms(formula)
    form_preds = terms(rlang::new_formula(lhs=NULL, rhs=forms$predictors))
    form_combined = rlang::new_formula(forms$outcome, expr(!!forms$predictors + !!forms$covariates))

    bp = hardhat::new_default_formula_blueprint(
        intercept = FALSE,
        composition = "matrix",
        ei_x = attr(form_preds, "term.labels"),
        ei_wgt = check_make_weights(!!enquo(weights), data),
        bounds = bounds,
        subclass = "ei_blueprint"
    )

    processed <- hardhat::mold(form_combined, data, blueprint=bp)
    ei_model_bridge(processed, ...)
}


#' @export
#' @rdname ei_model
ei_model.ei_spec <- function(x, bounds=NULL, ...) {
    spec = x
    x = spec[c(attr(spec, "ei_x"), attr(spec, "ei_z"))]
    y = spec[attr(spec, "ei_y")]
    bp = hardhat::new_default_xy_blueprint(
        intercept = FALSE,
        composition = "matrix",
        ei_x = attr(spec, "ei_x"),
        ei_wgt = attr(spec, "ei_wgt"),
        bounds = bounds,
        subclass = "ei_blueprint"
    )

    processed = hardhat::mold(x, y, blueprint=bp)
    ei_model_bridge(processed, ...)
}


#' @export
#' @rdname ei_model
ei_model.data.frame <- function(x, y, z, weights, bounds=NULL, ...) {
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
        ei_wgt = check_make_weights(weights, NULL),
        bounds = bounds,
        subclass = "ei_blueprint"
    )
    x = cbind(x, z)

    processed <- hardhat::mold(x, y, blueprint=bp)
    ei_model_bridge(processed, ...)
}

#' @export
#' @rdname ei_model
ei_model.matrix <- function(x, y, z, weights, bounds=NULL, ...) {
    ei_model.data.frame(x, y, z, weights, bounds=bounds, ...)
}


#' @export
#' @rdname ei_model
ei_model.default <- function(x, ...) {
    if (missing(x))
        cli_abort("{.fn ei_model} requires arguments.", call=NULL)
    cli_abort("{.fn ei_model} is not defined for a {.cls {class(x)}}.", call=NULL)
}

# Creates bounds _after_ molding outcomes
#' @exportS3Method hardhat::run_mold
run_mold.ei_blueprint <- function(blueprint, ...) {
    processed = NextMethod("run_mold")

    # update bounds
    bounds = ei_bounds(processed$blueprint$bounds, processed$outcomes)
    processed$blueprint = hardhat::update_blueprint(
        processed$blueprint,
        bounds = bounds
    )

    processed
}


# Bridge ----------------------------------------------------------------------

ei_model_bridge <- function(processed, ...) {
  x = processed$predictors
  idx_x = match(processed$blueprint$ei_x, colnames(x))
  z = x[, -idx_x, drop=FALSE]
  x = x[, idx_x, drop=FALSE]
  if (ncol(x) == 1) {
      x = cbind(x, 1 - x)
      colnames(x)[2] = ".other"
  }
  y = as.matrix(processed$outcomes)

  fit <- ei_model_impl(x, y, z, processed$blueprint$ei_wgt, processed$blueprint$bounds)

  do.call(new_ei_model, c(fit, list(blueprint = processed$blueprint)))
  # new_ei_model(
  #   coefs = fit$coefs,
  #   blueprint = processed$blueprint
  # )
}

# Model type ------------------------------------------------------------------

new_ei_model <- function(..., blueprint) {
    hardhat::new_model(..., blueprint = blueprint, class = "ei_model")
}

#' @export
print.ei_model <- function(x, ...) {
    cat_line("An ecological inference model")
    cat_line("<placeholder print method>")
    if ("b_global" %in% names(x)) {
        cat_line("Global estimate:")
        print(x$b_global)
    } else {
        cat_line("Eta estimate:")
        print(x$est$eta)
    }
}
