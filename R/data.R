#' 1968 U.S. Presidential Election Data
#'
#' County-level dataset containing election results and demographics from the
#' 1968 U.S. presidential election in the states of Virginia, North Carolina,
#' South Carolina, Georgia, Florida, Alabama, Mississippi, Louisiana, Arkansas,
#' Tennessee, and Texas.
#'
#' @format A data frame with 1,143 rows and 41 variables:
#' \describe{
#'   \item{`fips`}{County FIPS code}
#'   \item{`state`}{State name}
#'   \item{`abbr`}{State abbreviation}
#'   \item{`region`}{Census region}
#'   \item{`division`}{Census division}
#'   \item{`county`}{County name}
#'   \item{`pop`}{Total population at the 1970 census}
#'   \item{`pop_city`, `pop_urban`, `pop_rural`}{Proportion of the population in city, urban and rural areas}
#'   \item{`pop_white`, `pop_black`, `pop_aian`, `pop_asian`, `pop_hisp`}{
#'     Proportion of the population in each racial group.
#'     The first four columns sum to 1, but `pop_hisp` should be considered
#'     separately.}
#'   \item{`vap`}{Voting-age population at the 1970 census}
#'   \item{`vap_white`, `vap_black`, `vap_other`}{
#'     Proportion of the voting-age population in each racial group}
#'   \item{`farm`, `nonfarm`}{Proportion of the population in farm and non-farm households}
#'   \item{`educ_elem`, `educ_hsch`, `educ_coll`}{Proportion of the population with
#'     elementary, high school, and college education}
#'   \item{`cvap`}{Citizen voting-age population at the 1970 census}
#'   \item{`cvap_white`, `cvap_black`, `cvap_other`}{
#'     Proportion of the citizen voting-age population in each racial group}
#'   \item{`inc_00_03k`, `inc_03_08k`, `inc_08_25k`, `inc_25_99k`}{
#'     Proportion of the population in households with income in each bracket.
#'     Median household income in 1970 was $9,870}
#'   \item{`pres_dem_hum`, `pres_rep_nix`, `pres_ind_wal`, `pres_abs`, `pres_oth`}{
#'     Proportion of votes for Humphrey, Nixon, and Wallace, and other candidates.
#'     `pres_abs` is calculated as one minus the Humphrey, Nixon, and Wallace vote shares.}
#'   \item{`pres_total`}{Total number of votes cast for president.}
#'   \item{`pres_turn`}{Proportion of the voting-age population that turned out to vote.}
#'   \item{`pres_turn_icpsr`}{Proportion of the voting-age population that
#'     turned out to vote using ICPSR data in the denominator.}
#'   \item{`adj`}{A FIPS-indexed adjacency list capturing the geographic (rook)
#'     adjacency between counties. To convert to a 1-indexed adjacency list,
#'     run `lapply(elec_1968$adj, function(x) match(x, elec_1968$fips))`.}
#' }
#'
#' @source
#' **Census data**:
#' Steven Manson, Jonathan Schroeder, David Van Riper, Katherine Knowles, Tracy
#' Kugler, Finn Roberts, and Steven Ruggles. IPUMS National Historical
#' Geographic Information System: Version 19.0. 1970 Decennial Census.
#' Minneapolis, MN: IPUMS. 2024. <http://doi.org/10.18128/D050.V19.0>
#'
#' **Presidential election data**:
#' Clubb, Jerome M., Flanigan, William H., and Zingale, Nancy H. Electoral Data
#' for Counties in the United States: Presidential and Congressional Races,
#' 1840-1972. Inter-university Consortium for Political and Social Research
#' (distributor), 2006-11-13. <https://doi.org/10.3886/ICPSR08611.v1>
#'
#'
"elec_1968"
