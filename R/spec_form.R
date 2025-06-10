#' Specify an ecological inference problem
#'
#' Uses tidy-select syntax to specify outcomes, predictors, and covariates.
#' The result of this function can be passed directly into [ei_ridge()] or
#' [ei_riesz()], or plotted with [`plot()`][plot.ei_spec].
#'
#' The function is lightweight and does not perform any checking of the
#' arguments, bounds, sum constraints, etc.  All of these checks are performed
#' by [ei_ridge()] or [ei_riesz()].
#'
#' @param data A data frame.
#' @param predictors <[`tidy-select`][dplyr::select]> Predictor variables.
#'   This is the `x` variable in ecological regression that is of primary interest.
#'   For example, the columns containing the percentage of each racial group.
#' @param outcome <[`tidy-select`][dplyr::select]> Outcome variables.
#'   This is the `y` variable in ecological regression that is of primary interest.
#'   For example, the columns containing the percentage of votes for each party.
#' @param total <[`tidy-select`][dplyr::select]> A variable containing the total
#'   number of observations in each aggregate unit. For example, the column
#'   containing the total number of voters. Required by default.
#' @param covariates <[`tidy-select`][dplyr::select]> Covariates.
#' @param strip Whether to strip common prefixes from column names within each group.
#'   For example, columns named `vap_white`, `vap_black`, and `vap_hisp` would be
#'   renamed `white`, `black` and `other` in the model and output.
#'
#' @returns An `ei_spec` object, which is a data frame with additional
#'   attributes recording `predictors`, `outcomes`, `total`, and `covariates`.
#'
#' @examples
#' data(elec_1968)
#' ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth, pres_total)
#'
#' @export
ei_spec = function(data, predictors, outcome, total, covariates=NULL, strip=TRUE) {
    predictors = try_fetch(
        eval_select(enquo(predictors), data, allow_empty=FALSE),
        error = function(cnd) rlang::abort("Predictor specification failed.", parent=cnd)
    )
    outcome =  try_fetch(
        eval_select(enquo(outcome), data, allow_empty=FALSE),
        error = function(cnd) rlang::abort("Outcome specification failed.", parent=cnd)
    )
    covariates = try_fetch(
        eval_select(enquo(covariates), data),
        error = function(cnd) rlang::abort("Covariate specification failed.", parent=cnd)
    )
    total = check_make_weights(!!enquo(total), data)


    if (isTRUE(strip)) {
        names(predictors) = str_strip_prefix(names(predictors))
        names(outcome) = str_strip_prefix(names(outcome))
        names(covariates) = str_strip_prefix(names(covariates))
    }
    cols = c(predictors, outcome, covariates)

    new_tibble(
        setNames(data[cols], names(cols)),
        ei_x = names(predictors),
        ei_y = names(outcome),
        ei_z = names(covariates),
        ei_n = total,
        class="ei_spec"
    )
}

#' @export
print.ei_spec = function(x, ..., n=5) {
    cat_line("EI Specification")
    covs = attr(x, "ei_z")
    covs_lbl = if (length(covs) == 0) "none" else "{.var {covs}}"
    cli::cat_bullet(c(
        format_inline("Predictors: {.var {attr(x, 'ei_x')}}"),
        format_inline("Outcome: {.var {attr(x, 'ei_y')}}"),
        format_inline(paste0("Covariates: ", covs_lbl))
    ))
    NextMethod(n=n)
}

#' @describeIn ei_spec Extract the totals from a specification
#' @param object An [ei_spec] object.
#' @param normalize If `TRUE`, normalize the totals to have mean 1.
#' @param ... Additional arguments (ignored).
#' @export
weights.ei_spec <- function(object, normalize = TRUE, ...) {
    n = attr(object, "ei_n")
    if (is.null(n)) {
        n = rep(1, nrow(object))
    }
    if (isTRUE(normalize)) {
        n / mean(n)
    } else{
        n
    }
}

#' Plot an EI specification
#'
#' Useful for quickly visualizing scatterplots of outcome versus predictor
#' variables.
#'
#' @param x An [ei_spec] object.
#' @param ... Additional arguments passed to [pairs()].
#' @param pch,cex As in [plot()]
#'
#' @examples
#' data(elec_1968)
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth, pres_total)
#' plot(spec)
#'
#' @export
plot.ei_spec = function(x, ..., pch=16, cex=0.2) {
    nm_x = attr(x, "ei_x")
    nm_y = attr(x, "ei_y")
    size = sqrt(attr(x, "ei_n"))
    size = size / mean(size)

    graphics::pairs(x, horInd=match(nm_y, names(x)), verInd=match(nm_x, names(x)),
                    pch=pch, cex=cex*size, ...)
    adj_y = seq(0, 1, 1/(length(nm_y)*2))[2*seq_along(nm_y)]
    adj_x = seq(0, 1, 1/(length(nm_x)*2))[2*seq_along(nm_x)]
    graphics::mtext(rev(nm_y), side=2, line=3, adj=adj_y, padj=1, font=2)
    graphics::mtext(nm_x, side=3, line=3, adj=adj_x, padj=1, font=2)
}


# Helper to parse EI formulas
ei_forms = function(formula) {
    rhs = f_rhs(formula)
    if (rlang::is_symbol(rhs) || !as.character(rhs[[1]]) == "|") {
        rhs = rlang::expr(!!rhs | 1)
    }

    c(
        outcome = f_lhs(formula),
        predictors = f_lhs(rhs),
        covariates = f_rhs(rhs)
    )
}

# Helper to infer bounds
ei_bounds = function(bounds, outcome) {
    if (is.null(bounds)) {
        lt_1 = all(outcome <= 1)
        gt_0 = all(outcome >= 0)
        if (lt_1 && gt_0) {
            c(0, 1)
        } else if (gt_0) {
            c(0, Inf)
        } else {
            c(-Inf, Inf)
        }
    } else if (is.numeric(bounds)) {
        if (length(bounds) != 2)
            cli_abort("{.arg bounds} must be a length-2 vector.", call=parent.frame())
        if (bounds[1] >= bounds[2])
            cli_abort("Lower bound must be strictly less than upper bound.", call=parent.frame())

        if (any(outcome < bounds[1]) || any(outcome > bounds[2])) {
            cli_abort("Some outcomes are outside the specified bounds of
                      [{bounds[1]}, {bounds[2]}].", call=parent.frame())
        }
        bounds
    } else {
        cli_abort("Improper {.arg bounds} specification: {.val {bounds}}", call=parent.frame())
    }
}

# Shared helper for checking weights/total argument
check_make_weights = function(x, data, n, arg = "total", required = TRUE) {
    if (rlang::quo_is_missing(enquo(x))) {
        if (isTRUE(required)) {
            cli_abort(c(
                "The {.arg {arg}} argument is required.",
                "i"="{.arg {arg}} should contain the total number of individuals in each unit.",
                ">"="To use uniform values (not recommended), pass {.arg {arg} = FALSE}."
            ), call=parent.frame())
        } else {
            x = FALSE
        }
    }
    x = eval_tidy(enquo(x), data)

    if (!is.null(data)) {
        n = nrow(data)
    }

    if (isFALSE(x)) {
        x = rep(1, n)
    }
    if (any(x < 0)) {
        cli_abort("Negative totals and weights are not allowed.", call=parent.frame())
    }

    if (length(x) != n) {
        cli_abort("Length of {.arg {arg}} ({length(x)}) does not match
                  the number of data observations ({n}).", call=parent.frame())
    }

    x
}


# Helper to remove common prefixes from strings
str_strip_prefix = function(x) {
    if (length(x) <= 1) return(x)
    if (all(nchar(x) == 0)) return(x)
    n = max(nchar(x), na.rm=TRUE)
    found = FALSE
    for (i in seq_len(n)) {
        if (length(unique(substr(x, i, i))) > 1) {
            found = TRUE
            break
        }
    }
    if (!found) {
        rep(NA_character_, length(x))
    } else {
        substr(x, i, n)
    }
}
