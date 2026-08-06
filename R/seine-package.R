#' @keywords internal
#'
#' @references
#' McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
#' estimation of conditional means from aggregate data.
#' Working paper [arXiv:2509.20194](https://arxiv.org/abs/2509.20194).
#'
#' Kuriwaki, S., & McCartan, C. (2026). The role of confounders and linearity
#' in ecological inference: A reassessment.
#' Working paper [arXiv:2601.07668](https://arxiv.org/abs/2601.07668).
"_PACKAGE"

## usethis namespace: start
#' @rawNamespace import(stats, except = c("filter", "lag"))
#' @importFrom cli cat_line cli_abort cli_warn cli_inform format_inline
#' @importFrom rlang check_dots_empty0
#' @importFrom rlang eval_tidy expr enquo f_lhs f_rhs `f_lhs<-` try_fetch
#' @importFrom tibble as_tibble new_tibble
#' @importFrom tidyselect eval_select
#' @useDynLib seine, .registration = TRUE
## usethis namespace: end
NULL
