#' Compute bounds on local ecological quantities
#'
#' For each observation, computes the minimum and maximum value of each local
#' estimand that is consistent with the accounting constraint \eqn{B %*% x = y},
#' where \eqn{B} is the matrix of local estimands, \eqn{x} are the predictors,
#' and \eqn{y} are the outcomes.
#'
#' @inheritParams ei_ridge
#' @inheritParams ei_est_local
#' @param formula A formula such as `y ~ x0 + x1` specifying the outcome `y`
#'   and the predictors of interest `x`. The predictors should form a partition,
#'   that is, `x0 + x1 = 1` for each observation.
#'   Users can be include more than two predictors as well, e.g.
#'   `pct_white + pct_black + pct_hisp + pct_other`.
#'   If there are just two predictors, it is acceptable to only include one in
#'   the formula; the other will be formed as 1 minus the provided predictor.
#' @param total <[`data-masking`][rlang::args_data_masking]> A variable
#'   containing the total number of observations in each aggregate unit. For
#'   example, the column containing the total number of voters. Required for
#'   computing weights unless `x` is an [ei_spec()] object.
#' @param global If `TRUE`, aggregate the bounds across units to produce bounds
#'   on the global estimands.
#' @param ... Not currently used, but required for extensibility.
#'
#' @returns A data frame with bounds. The `.row` column in the output
#'   corresponds to the observation index in the input. The `min` and `max`
#'   columns contain the minimum and maximum values for each local estimand.
#'   The `wt` column contains the product of the predictor variable and total
#'   for each observation. Taking a weighted average of the bounds against this
#'   column will produce global bounds. It has class `ei_bounds`.
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
#'                total = pres_total, covariates = c(state, pop_urban, farm))
#'
#' ei_bounds(spec, bounds = c(0, 1))
#' ei_bounds(spec, bounds = c(0, 1), global = TRUE)
#'
#' # Infer bounds
#' ei_bounds(pres_ind_wal ~ vap_white, data = elec_1968, total = pres_total, bounds = NULL)
#'
#' # manually aggregate min/max
#' # easier with dplyr:
#' # summarize(across(min:max, ~ weighted.mean(.x, wt)), .by=c(predictor, outcome))
#' grp_units = split(ei_bounds(spec, bounds = c(0, 1)), ~ predictor + outcome)
#' do.call(rbind, lapply(grp_units, function(b) {
#'     tibble::tibble(
#'         predictor = b$predictor[1],
#'         outcome = b$outcome[1],
#'         min = weighted.mean(b$min, b$wt),
#'         max = weighted.mean(b$max, b$wt)
#'     )
#' }))
#'
#' @export
ei_bounds <- function(x, ..., total, contrast = NULL, bounds = c(0, 1), global = FALSE) {
    UseMethod("ei_bounds")
}

#' @export
#' @rdname ei_bounds
ei_bounds.ei_spec <- function(x, total, contrast = NULL, bounds = c(0, 1), global = FALSE, ...) {
    spec = x
    validate_ei_spec(spec)

    x_nm = attr(spec, "ei_x")
    x_mat = as.matrix(pull_x(spec, x_nm))
    check_preds(x_mat, call=parent.frame())
    y_nm = attr(spec, "ei_y")
    y_mat = as.matrix(spec[, y_nm, drop = FALSE])

    bounds = check_bounds(bounds, y_mat)

    if (missing(total)) {
        total = as.numeric(attr(spec, "ei_n"))
    } else {
        total = as.numeric(eval_tidy(enquo(total), spec))
    }

    ei_bounds_bridge(x_mat, y_mat, total, contrast, bounds, global)
}

#' @export
#' @rdname ei_bounds
ei_bounds.formula <- function(formula, data, total, contrast = NULL, bounds = c(0, 1), global = FALSE, ...) {
    forms = ei_forms(formula)
    form_preds = terms(rlang::new_formula(lhs = NULL, rhs = forms$predictors))
    form_out = terms(rlang::new_formula(forms$outcome, rhs = NULL))

    bp = hardhat::new_default_formula_blueprint(
        intercept = FALSE,
        composition = "matrix",
        indicators = "one_hot",
        bounds = bounds,
        subclass = "ei_ridge_blueprint"
    )

    processed <- hardhat::mold(form_out, data, blueprint=bp)

    x_nm = attr(form_preds, "term.labels")
    x = as.matrix(pull_x(data, x_nm))
    check_preds(x, call=parent.frame())

    y = as.matrix(processed$outcomes)

    total = as.numeric(eval_tidy(enquo(total), data))

    ei_bounds_bridge(x, y, total, contrast, processed$blueprint$bounds, global)
}


#' @export
#' @rdname ei_bounds
ei_bounds.data.frame <- function(x, y, total, contrast = NULL, bounds = c(0, 1), global = FALSE, ...) {
    x_mat = as.matrix(x)
    check_preds(x_mat, call = rlang::new_call(rlang::sym("ei_bounds")))
    y_mat = as.matrix(y)

    # Get total
    if (missing(total)) {
        cli_abort("{.arg total} is required when using data frames or matrices.")
    }
    total = as.numeric(total)

    bounds = check_bounds(bounds, y_mat)

    ei_bounds_bridge(x_mat, y_mat, total, contrast, bounds, global)
}

#' @export
#' @rdname ei_bounds
ei_bounds.matrix <- function(x, y, total, contrast = NULL, bounds = c(0, 1), global = FALSE, ...) {
    ei_bounds.data.frame(x, y, total, contrast, bounds, global, ...)
}

#' @export
#' @rdname ei_bounds
ei_bounds.default <- function(x, ...) {
    if (missing(x)) {
        cli_abort("{.fn ei_bounds} requires arguments.", call = NULL)
    }
    cli_abort("{.fn ei_bounds} is not defined for a {.cls {class(x)}}.", call = NULL)
}

# Implementation --------------------------------------------------------------

ei_bounds_bridge <- function(x, y, total, contrast, bounds, global = FALSE) {
    n = nrow(x)
    n_x = ncol(x)
    n_y = ncol(y)

    if (identical(bounds, c(-Inf, Inf))) {
        cli_abort("At least one bound must be provided for {.fn ei_bounds}.", call=parent.frame())
    }
    if (any(is.na(x))) cli_abort("Missing values found in predictors.", call=parent.frame())
    if (any(is.na(y))) cli_abort("Missing values found in outcome.", call=parent.frame())

    result = ei_bounds_impl(x, y, total, contrast, bounds)

    x_nm = colnames(x)
    y_nm = colnames(y)

    if (isTRUE(global)) {
        wt = total * (matrix(1, 1, n_y) %x% x)
        wt = scale_cols(wt, 1 / colSums(wt))
        tibble::new_tibble(
            list(
                predictor = rep(x_nm, n_y),
                outcome = rep(y_nm, each = n_x),
                min = colSums(result$min * wt),
                max = colSums(result$max * wt)
            ),
            class = "ei_bounds"
        )
    } else {
        tibble::new_tibble(
            list(
                .row = rep(seq_len(n), n_x * n_y),
                predictor = rep(rep(x_nm, each = n), n_y),
                outcome = rep(y_nm, each = n * n_x),
                wt = if (is.null(contrast)) rep(c(x * total), n_y) else NULL,
                min = c(result$min),
                max = c(result$max)
            ),
            class = "ei_bounds"
        )
    }
}

ei_bounds_impl <- function(x, y, total, contrast, bounds) {
    if (!is.null(contrast)) {
        cli_abort("{.arg contrast} is not yet implemented for {.fn ei_bounds}.")
    }

    R_bounds_lp(x, y, as.double(bounds))
}

#' @describeIn ei_bounds Format bounds as an array with dimensions
#'   `<rows>*<predictors>*<outcomes>*2`. Does not work if the object has been sorted.
#' @param x An object of class `ei_bounds`
#' @param ... Additional arguments (ignored)
#' @export
as.array.ei_bounds = function(x, ...) {
    nm_x = unique(x$predictor)
    nm_y = unique(x$outcome)
    n_x = length(nm_x)
    n_y = length(nm_y)
    n = nrow(x) / n_x / n_y

    out = array(
        c(x$min, x$max),
        dim = c(n, n_x, n_y, 2),
        dimnames = list(NULL, nm_x, nm_y, c("min", "max"))
    )
    out
}
