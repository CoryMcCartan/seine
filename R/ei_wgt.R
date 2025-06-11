#' Estimation weighting models
#'
#' Several built-in helper functions to generate estimation weights from a
#' vector of unit totals, or an existing [ei_spec()] object.
#'
#' @param x A numeric vector of unit totals, or an existing [ei_spec()] object.
#'
#' @returns A numeric vector of estimation weights with the same number of
#'   observations as `x`. These will have mean 1.
#'
#' @examples
#' data(elec_1968)
#'
#' ei_wgt_unif(elec_1968$pres_total)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal, total=pres_total)
#' ei_wgt_prop(spec)
#' ei_wgt_sqrt(spec)
#'
#' @name ei_wgt
NULL

#' @describeIn ei_wgt Uniform weights. Appropriate if the unit-level variance
#'   is constant, i.e., homosekdastic.
#' @export
ei_wgt_unif <- function(x) {
    x = get_tot(x)
    rep(1, length(x))
}

#' @describeIn ei_wgt Weights proportional to the totals. Appropriate if the
#'   unit-level variance is inversely proportional to the number of
#'   observations.
#' @export
ei_wgt_prop <- function(x) {
    x = get_tot(x)
    x / mean(x)
}

#' @describeIn ei_wgt Weights proportional to the square root of the totals.
#'   Appropriate if the unit-level variance is inversely proportional to the
#'   square root of the number of observations.
#' @export
ei_wgt_sqrt <- function(x) {
    x = sqrt(get_tot(x))
    x / mean(x)
}

get_tot <- function(x) {
    if (is.numeric(x)) {
        if (any(x < 0)) {
            cli_abort("{.arg x} must be nonnegative.", call = parent.frame())
        }
        x
    } else if (inherits(x, "ei_spec")) {
        attr(x, "ei_n")
    } else {
        cli_abort("{.arg x} must be a numeric vector or {.cls ei_spec} object.",
                  call = parent.frame())
    }
}
