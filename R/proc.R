# Data pre- and post-processing

#' Convert counts to proportions
#'
#' Divides counts in specified columns by a specified total or by their sum,
#' possibly storing the result in a new column.
#' Also checks for out-of-bounds and missing values, and can create a column
#' containing the remainder amount so that all proportions sum to 1.
#'
#' @param data A data frame.
#' @param ... <[`tidy-select`][tidyselect::select_helpers]> Columns to convert
#'   to proportions. If one column is completely missing, it will be imputed as
#'   1 minus the total.
#' @param .total <[`tidy-select`][tidyselect::select_helpers]> Column to use as
#'   the total. If this column does not exist, it will be created as the sum of
#'   the selected columns. Supports renaming syntax.
#' @param .other <[`tidy-select`][tidyselect::select_helpers]> Column to store
#'   the remainder, so that the selected columns plus `.other` sum to 1. If the
#'   selected columns do sum to 1 after normalization, this argument will not be
#'   used; otherwise it will be created or overwritten. The calculation of
#'   `.other` is performed _after_ clamping (see below).
#' @param clamp Proportions that are `clamp` below 0 or above 1 will be rounded
#'   to 0 and 1, respectively.  Values outside `clamp` will throw an error.
#'   Set `clamp=0` to disable or `clamp=Inf` to allow for out-of-bounds
#'   proportions (not recommended).
#'
#' @returns A modified data frame.  Unselected columns are unmodified.
#'
#' @examples
#' data(elec_1968)
#' # Make a data frame with counts
#' d_unnorm = with(head(elec_1968, 10), data.frame(
#'   vap = vap,
#'   vap_white = vap * vap_white,
#'   vap_black = vap * vap_black,
#'   vap_other = vap * vap_other
#' ))
#'
#' ei_proportions(d_unnorm, vap_white:vap_black, .total=vap) # `.other` column created
#' ei_proportions(d_unnorm, vap_white:vap_other) # no total provided
#' # renaming allowed
#' ei_proportions(d_unnorm, white=vap_white, black=vap_black,
#'                .total=c(total=vap), .other="vap_other")
#' @export
ei_proportions = function(data, ..., .total=".total", .other=".other", clamp=1e-3) {
    cols = eval_select(expr(c(...)), data, allow_empty=FALSE)
    .total = enquo(.total)
    .other = enquo(.other)
    # will error if column doesn't exist but handles character values and NULL OK
    total = if (rlang::quo_is_symbolic(.total)) eval_select(.total, data) else character(0)
    other = if (rlang::quo_is_symbolic(.other)) eval_select(.other, data) else character(0)

    # calculate totals
    if (length(total) == 0) {
        total = rowSums(data[cols])
        if (!rlang::quo_is_null(.total)) { # if not NULL, create column
            data[[rlang::quo_get_expr(.total)]] = total
        }
    } else if (length(total) == 1) {
        names(data)[total] = names(total)
        total = data[[total]]
    } else {
        cli_abort("{.arg total} must refer to a single column.")
    }

    data[cols] = data[cols] / total
    old_names = names(data)[cols]
    names(data)[cols] = names(cols)

    if (any(is.na(data[cols]))) {
        has_na = old_names[which(colSums(is.na(data[cols])) > 0)]
        cli_abort("Some selected columns have missing values: {.var {has_na}}.")
    }

    # clamp
    x = data[cols]
    x[-clamp < x & x < 0] = 0
    x[1 < x & x < 1 + clamp] = 1
    totals = rowSums(x)
    exc = 1 < totals & totals < 1 + clamp
    x[exc, ] = x[exc, ] / totals[exc]
    totals[exc] = 1
    data[cols] = x

    if (any(totals > 1)) {
        cli_abort(c("Selected columns total more than 1.",
                    "i"="Proportions should be nonnegative and sum to 1."))
    }
    if (any(x < -clamp)) {
        has_neg = old_names[which(colSums(x < -clamp) > 0)]
        cli_abort(c("Some selected columns have negative values: {.var {has_neg}}.",
                    "i"="Proportions should be nonnegative and sum to 1."))
    }

    # handle .other
    if (!isTRUE(all.equal(totals, rep(1, nrow(x))))) {
        leftover = pmax(1 - totals, 0)

        if (length(other) == 0) {
            if (rlang::quo_is_null(.other)) {
                cli_abort(c("{.arg .other} cannot be {.val NULL}.",
                            "i"="Selected columns do not sum to 1"))
            }
            data[[rlang::quo_get_expr(.other)]] = leftover
        } else {
            names(data)[other] = names(other)
            data[[other]] = leftover
        }
    }

    data
}
