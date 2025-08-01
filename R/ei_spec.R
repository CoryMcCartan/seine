#' Specify an ecological inference problem
#'
#' Uses tidy-select syntax to specify outcomes, predictors, and covariates.
#' The result of this function can be passed directly into [ei_ridge()]
#' or [ei_riesz()], or plotted with [`plot()`][plot.ei_spec].
#'
#' The function is lightweight and does not perform any checking of the
#' arguments, bounds, sum constraints, etc.  All of these checks are performed
#' by functions that use `ei_spec` objects.
#'
#' @param data A data frame.
#' @param predictors <[`tidy-select`][tidyselect::select_helpers]> Predictor
#'   variables. This is the `x` variable in ecological regression that is of
#'   primary interest. For example, the columns containing the percentage of
#'   each racial group.
#' @param outcome <[`tidy-select`][tidyselect::select_helpers]> Outcome
#'   variables. This is the `y` variable in ecological regression that is of
#'   primary interest. For example, the columns containing the percentage of
#'   votes for each party.
#' @param total <[`tidy-select`][tidyselect::select_helpers]> A variable
#'   containing the total number of observations in each aggregate unit. For
#'   example, the column containing the total number of voters. Required by
#'   default.
#' @param covariates <[`tidy-select`][tidyselect::select_helpers]> Covariates.
#' @param strip Whether to strip common prefixes from column names within each group.
#'   For example, columns named `vap_white`, `vap_black`, and `vap_hisp` would be
#'   renamed `white`, `black` and `other` in the model and output.
#'
#' @returns An `ei_spec` object, which is a data frame with additional
#'   attributes recording `predictors`, `outcomes`, `total`, and `covariates`.
#'
#' @examples
#' data(elec_1968)
#' ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs, pres_total)
#'
#' @export
ei_spec = function(data, predictors, outcome, total, covariates=NULL, strip=FALSE) {
    predictors = try_fetch(
        eval_select(enquo(predictors), data, allow_empty=FALSE),
        error = function(cnd) cli_abort("Predictor specification failed.", parent=cnd)
    )
    outcome =  try_fetch(
        eval_select(enquo(outcome), data, allow_empty=FALSE),
        error = function(cnd) cli_abort("Outcome specification failed.", parent=cnd)
    )
    covariates = try_fetch(
        eval_select(enquo(covariates), data),
        error = function(cnd) cli_abort("Covariate specification failed.", parent=cnd)
    )
    total = check_make_weights(!!enquo(total), data)


    if (isTRUE(strip)) {
        names(predictors) = str_strip_prefix(names(predictors))
        names(outcome) = str_strip_prefix(names(outcome))
        names(covariates) = str_strip_prefix(names(covariates))
    }
    cols = c(predictors, outcome, covariates)

    new_ei_spec(
        data = setNames(data[cols], names(cols)),
        x = names(predictors),
        y = names(outcome),
        z = names(covariates),
        n = total
    )
}

# Internal constructor
new_ei_spec = function(data, x, y, z, n, ...) {
    new_tibble(data, ei_x = x, ei_y = y, ei_z = z, ei_n = n, class="ei_spec")
}

# Internal validator
validate_ei_spec = function(x) {
    if (!inherits(x, "ei_spec")) {
        cli_abort("Object must be an {.cls ei_spec} object.", call=parent.frame())
    }
    if (!is.data.frame(x)) {
        cli_abort("An {.cls ei_spec} object must be a data frame.", call=parent.frame())
    }
    if (is.null(attr(x, "ei_x"))) {
        cli_abort("No predictors specified in {.cls ei_spec} object.", call=parent.frame())
    }
    if (is.null(attr(x, "ei_y"))) {
        cli_abort("No outcome specified in {.cls ei_spec} object.", call=parent.frame())
    }
    if (is.null(attr(x, "ei_n"))) {
        cli_abort("No total specified in {.cls ei_spec} object.", call=parent.frame())
    }

    if (!all(attr(x, "ei_x") %in% names(x))) {
        cli_abort(c(
            "Some predictors are not present in the data frame.",
            ">" = paste0("Missing: {.var ", setdiff(attr(x, "ei_x"), names(x)), "}")
        ), call=parent.frame())
    }
    if (!all(attr(x, "ei_y") %in% names(x))) {
        cli_abort(c(
            "Some outcomes are not present in the data frame.",
            ">" = paste0("Missing: {.var ", setdiff(attr(x, "ei_y"), names(x)), "}")
        ), call=parent.frame())
    }
    if (!all(attr(x, "ei_z") %in% names(x))) {
        cli_abort(c(
            "Some covariates are not present in the data frame.",
            ">" = paste0("Missing: {.var ", setdiff(attr(x, "ei_z"), names(x)), "}")
        ), call=parent.frame())
    }
    if (length(attr(x, "ei_n")) != nrow(x)) {
        cli_abort("Length of {.var ei_n} ({length(attr(x, 'ei_n'))}) does not
                  match the number of data observations ({nrow(x)}).",
                  call=parent.frame())
    }

    xm = as.matrix(x[, attr(x, "ei_x")])
    if (storage.mode(xm) != "double") {
        cli_abort("Predictors must be numeric; found {.cls {storage.mode(xm)}}.",
                  call=parent.frame())
    }
    check_preds(xm)

    invisible(x)
}

reconstruct_ei_spec = function(data, old) {
    if (!missing(old)) {
        if (is.null(attr(data, "ei_x"))) {
            attr(data, "ei_x") = attr(old, "ei_x")
        }
        if (is.null(attr(data, "ei_y"))) {
            attr(data, "ei_y") = attr(old, "ei_y")
        }
        if (is.null(attr(data, "ei_z"))) {
            attr(data, "ei_z") = attr(old, "ei_z")
        }
        if (is.null(attr(data, "ei_n"))) {
            attr(data, "ei_n") = attr(old, "ei_n")
        }
    }

    attr(data, "ei_z") = intersect(attr(data, "ei_z"), names(data))

    if (!inherits(data, "ei_spec")) {
        class(data) = c("ei_spec", class(data))
    }

    data
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
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs, pres_total)
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
ei_bounds = function(bounds, outcome, clamp=1e-3) {
    if (isFALSE(bounds)) return(c(-Inf, Inf))
    if (is.null(bounds)) {
        if (is.null(outcome)) outcome = c(-Inf, Inf)
        lt_1 = all(outcome <= 1 + clamp)
        gt_0 = all(outcome >= -clamp)
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

        if (any(outcome < bounds[1] - clamp) || any(outcome > bounds[2] + clamp)) {
            cli_warn("Some outcomes are outside the specified bounds of
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
    if (any(is.na(x))) {
        cli_abort("Missing values found in {.arg {arg}}.", call=parent.frame())
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

check_preds = function(x, tol = 1e-6, call=parent.frame()) {
    if (!isTRUE(all.equal(rowSums(x), rep(1, nrow(x)), tolerance = tol))) {
        cli_warn("Predictors should sum to 1 in every row.", call=call)
    }
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
