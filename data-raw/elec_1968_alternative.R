library(tidyverse)
library(haven)
library(sf)
library(rmapshaper)

# codes ---
stateicpsr_code <- c(
    Connecticut = 01,
    Maine = 02,
    Massachusetts = 03,
    `New Hampshire` = 04,
    `Rhode Island` = 05,
    Vermont = 06,
    Delaware = 11,
    `New Jersey` = 12,
    `New York` = 13,
    Pennsylvania = 14,
    Illinois = 21,
    Indiana = 22,
    Michigan = 23,
    Ohio = 24,
    Wisconsin = 25,
    Iowa = 31,
    Kansas = 32,
    Minnesota = 33,
    Missouri = 34,
    Nebraska = 35,
    `North Dakota` = 36,
    `South Dakota` = 37,
    Virginia = 40,
    Alabama = 41,
    Arkansas = 42,
    Florida = 43,
    Georgia = 44,
    Louisiana = 45,
    Mississippi = 46,
    `North Carolina` = 47,
    `South Carolina` = 48,
    Texas = 49,
    Kentucky = 51,
    Maryland = 52,
    Oklahoma = 53,
    Tennessee = 54,
    `West Virginia` = 56,
    Arizona = 61,
    Colorado = 62,
    Idaho = 63,
    Montana = 64,
    Nevada = 65,
    `New Mexico` = 66,
    Utah = 67,
    Wyoming = 68,
    California = 71,
    Oregon = 72,
    Washington = 73,
    Alaska = 81,
    Hawaii = 82) |>
    enframe(name = "state", value = "icpsr_state")

# https://www.nhgis.org/
cens1970 <- read_csv("data-raw/nhgis0003_csv/nhgis0003_ds94_1970_county.csv") |>
    janitor::clean_names()

samp_inc <- read_csv("data-raw/nhgis0003_csv/nhgis0003_ds99_1970_county.csv") |>
    janitor::clean_names()


shp1970 <- read_sf("data-raw/nhgis0004_shape/nhgis0004_shapefile_tl2008_us_county_1970/")

# https://www.icpsr.umich.edu/web/ICPSR/studies/8611
pres_county <- read_dta("data-raw/08611-0001_renamed.dta")


# format ICPSR election file ---
pres_out <- pres_county |>
    mutate(
        icpsr_state = as.numeric(V1),
    ) |>
    tidylog::left_join(stateicpsr_code) |>
    tidylog::left_join(ccesMRPprep::states_key, by = "state") |>
    select(state, st, starts_with("county"), region, division, matches("pre_1968"), icpsr_state) |>
    mutate(across(matches("_(hum|nix|wallace)"), \(x) x / 100)) |>
    tidylog::filter(!is.na(pre_1968_total)) |>
    mutate(countya = county_id/10, .after = county_id) # padded 0 for some reason

# Format NHGIS shape file ---
shp_abbrv <- shp1970 |>
    janitor::clean_names() |>
    # typos?
    tidylog::mutate(
        icpsrcty = replace(icpsrcty, icpsrnam == "JACKSON" & statenam == "Georgia", 1570),
        icpsrcty = replace(icpsrcty, icpsrnam == "NEOSHO/DORN" & statenam == "Kansas", 1330),
        icpsrcty = replace(icpsrcty, icpsrnam == "NESS" & statenam == "Kansas", 1350),
    ) |>
    transmute(
        icpsr_state = as.numeric(icpsrst),
        county_id = as.numeric(icpsrcty),
        nhgisst, nhgiscty, nhgisnam,
        geometry) |>
    st_simplify(dTolerance = 2000) # seems to be in meters?

# Format Census 1970 tabulation ---
cens_inc <- samp_inc |>
    mutate(total = rowSums(pick(starts_with("c3t")))) |>
    transmute(
        state,
        countya = as.numeric(countya),
        cen_1970_1k = c3t001/total,
        cen_1970_2k = rowSums(pick(matches("c3t00[1-2]")))/total,
        cen_1970_5k = rowSums(pick(matches("c3t00[1-5]")))/total,
        cen_1970_10k = rowSums(pick(matches("c3t0(0[1-9]|10)")))/total,
        cen_1970_25k = rowSums(pick(matches("c3t0(0[1-9]|10|11|12|13)")))/total,
    )

cens <- cens1970 |>
    transmute(state,
              county,
              countya = as.numeric(countya),
              areaname,
              cen_1970_total = cbc001,
              cen_1970_rural = ccn001 / cbc001,
              cen_1970_urban = cbg001 / cbc001,
              cen_1970_city = cbd001 / cbc001,
              cen_1970_white = cbw001 / cbc001,
              cen_1970_black = cbw002 / cbc001,
              cen_1970_indian = cbw003 / cbc001,
              cen_1970_raceother = (cbc001 - (cbw001 + cbw002 + cbw003))/cbc001) |>
    tidylog::left_join(cens_inc, by = c("state", "countya"))



# Join -----
county_dat <- pres_out |>
    tidylog::inner_join(
        cens,
        by = c("state", "countya"),
        relationship = "one-to-one"
    ) |>
    tidylog::inner_join(
        shp_abbrv,
        by = c("icpsr_state", "county_id"),
        relationship = "one-to-one"
    ) |>
    filter(!st %in% c("AK", "HI")) |>
    st_as_sf()



# Plot ---
ggplot(county_dat) +
    geom_sf(aes(fill = pre_1968_wallace), linewidth = 0) +
    wacolors::scale_fill_wa_c()
