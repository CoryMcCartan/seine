#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import stats
#' @importFrom cli cat_line cli_abort cli_warn cli_inform format_inline
#' @importFrom rlang check_dots_empty0
#' @importFrom rlang eval_tidy expr enquo f_lhs f_rhs `f_lhs<-` try_fetch
#' @importFrom tibble new_tibble
#' @importFrom tidyselect eval_select
#' @useDynLib seine, .registration = TRUE
## usethis namespace: end
NULL


rlang::on_load(rlang::local_use_cli())
