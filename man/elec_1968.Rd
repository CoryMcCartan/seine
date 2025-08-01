% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{elec_1968}
\alias{elec_1968}
\title{1968 U.S. presidential election data}
\format{
A data frame with 1,143 rows and 41 variables:
\describe{
\item{\code{fips}}{County FIPS code}
\item{\code{state}}{State name}
\item{\code{abbr}}{State abbreviation}
\item{\code{region}}{Census region}
\item{\code{division}}{Census division}
\item{\code{county}}{County name}
\item{\code{pop}}{Total population at the 1970 census}
\item{\code{pop_city}, \code{pop_urban}, \code{pop_rural}}{Proportion of the population in city, urban and rural areas}
\item{\code{pop_white}, \code{pop_black}, \code{pop_aian}, \code{pop_asian}, \code{pop_hisp}}{
Proportion of the population in each racial group.
The first four columns sum to 1, but \code{pop_hisp} should be considered
separately.}
\item{\code{vap}}{Voting-age population at the 1970 census}
\item{\code{vap_white}, \code{vap_black}, \code{vap_other}}{
Proportion of the voting-age population in each racial group}
\item{\code{farm}, \code{nonfarm}}{Proportion of the population in farm and non-farm households}
\item{\code{educ_elem}, \code{educ_hsch}, \code{educ_coll}}{Proportion of the population with
elementary, high school, and college education}
\item{\code{cvap}}{Citizen voting-age population at the 1970 census}
\item{\code{cvap_white}, \code{cvap_black}, \code{cvap_other}}{
Proportion of the citizen voting-age population in each racial group}
\item{\code{inc_00_03k}, \code{inc_03_08k}, \code{inc_08_25k}, \code{inc_25_99k}}{
Proportion of the population in households with income in each bracket.
Median household income in 1970 was $9,870}
\item{\code{pres_dem_hum}, \code{pres_rep_nix}, \code{pres_ind_wal}, \code{pres_abs}, \code{pres_oth}}{
Proportion of votes for Humphrey, Nixon, and Wallace, and other candidates.
\code{pres_abs} is calculated as one minus the Humphrey, Nixon, and Wallace vote shares.}
\item{\code{pres_total}}{Total number of votes cast for president.}
\item{\code{pres_turn}}{Proportion of the voting-age population that turned out to vote.}
\item{\code{pres_turn_icpsr}}{Proportion of the voting-age population that
turned out to vote using ICPSR data in the denominator.}
\item{\code{adj}}{A FIPS-indexed adjacency list capturing the geographic (rook)
adjacency between counties. To convert to a 1-indexed adjacency list,
run \code{lapply(elec_1968$adj, function(x) match(x, elec_1968$fips))}.}
}
}
\source{
\strong{Census data}:
Steven Manson, Jonathan Schroeder, David Van Riper, Katherine Knowles, Tracy
Kugler, Finn Roberts, and Steven Ruggles. IPUMS National Historical
Geographic Information System: Version 19.0. 1970 Decennial Census.
Minneapolis, MN: IPUMS. 2024. \doi{10.18128/D050.V19.0}

\strong{Presidential election data}:
Clubb, Jerome M., Flanigan, William H., and Zingale, Nancy H. Electoral Data
for Counties in the United States: Presidential and Congressional Races,
1840-1972. Inter-university Consortium for Political and Social Research
(distributor), 2006-11-13. \doi{10.3886/ICPSR08611.v1}
}
\usage{
elec_1968
}
\description{
County-level dataset containing election results and demographics from the
1968 U.S. presidential election in the states of Virginia, North Carolina,
South Carolina, Georgia, Florida, Alabama, Mississippi, Louisiana, Arkansas,
Tennessee, and Texas.
}
\keyword{datasets}
