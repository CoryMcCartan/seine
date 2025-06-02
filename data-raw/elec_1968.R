library(tidyverse)
library(sf)
#'  * 1970 Decennial Census: Census and shapefiles from <https://www.nhgis.org/>
#     At county level
#     - NT103, NT114, NT75, NT25 (breakdowns: all 4 race/eth categories)
#     - A00, AA1, AC5, A35, A38, A57, B18
#'  * Presidential election data from <https://www.icpsr.umich.edu/web/ICPSR/studies/8611>
#'  * Also see: Black and Black (1973), The Wallace Vote in Alabama,
#'     JOP, <http://www.jstor.com/stable/2129154>;
#'     Wright (1977), Contextual Models of Electoral Behavior:
#'     The Southern Wallace Vote, APSR, <https://www.jstor.org/stable/1978344>
# Will need to update file names depending on extract number

c1970_raw1 <- read_csv("data-raw/nhgis0004_csv/nhgis0004_ds98_1970_county.csv") |>
  janitor::clean_names()
c1970_raw2 <- read_csv("data-raw/nhgis0004_csv/nhgis0004_ds99_1970_county.csv") |>
  janitor::clean_names()
c1970_raw3 <- read_csv("data-raw/nhgis0004_csv/nhgis0004_ts_nominal_county.csv") |>
  janitor::clean_names()

c1970_raw1 <- c1970_raw1 |>
  select(statea, countya, farm = c0v001, nonfarm = c0v002, c06001:c06010) |>
  mutate(
    educ_elem = rowSums(pick(c06001:c06005)),
    educ_hsch = rowSums(pick(c06006:c06007)),
    educ_coll = rowSums(pick(c06008:c06010))
  ) |>
  select(-c06001:-c06010)

c1970_raw2 <- c1970_raw2 |>
  transmute(
    statea = statea, countya = countya,
    # these CVAP are based on 21+ and so will be a bit high
    cvap = c12005 + c12010,
    cvap_white = na_if(c12aa005, -1) + na_if(c12aa010, -1),
    cvap_black = na_if(c12ab005, -1) + na_if(c12ab010, -1),
    # cvap_hisp = na_if(c12ac005, -1) + na_if(c12ac010, -1),
    cvap_denom = cvap + c12015, # an estimate of VAP built as `cvap` is
    # median inc in 1970 was around $8k, poverty line was $3100 for a family
    inc_00_03k = rowSums(pick(matches("c3t00[1-3]"))),
    inc_03_08k = rowSums(pick(matches("c3t00[4-8]"))),
    inc_08_25k = rowSums(pick(matches("c3t0(08|09|10|11|12|13)"))),
    inc_25_99k = rowSums(pick(matches("c3t01[4-5]")))
  )

# In general: voting age in 1968 was 21 except in Georgia and Kentucky (18) and
# Alaska and Hawaii (but excluded here)
# See <https://www2.census.gov/library/publications/1968/demographics/P25-406.pdf>
#
# So in April 1970 you had to be ~ 23 (most of US) or 20 (GA + KY) to have voted in 1968
c1970_raw3 <- c1970_raw3 |>
  mutate(across(c(a38au1970:a38dw1970, aa1au1970:ac5dw1970), ~ if_else(. >= 0, ., NA_real_))) |>
  transmute(
    state = state, county = county, statea = statefp, countya = countyfp,
    pop = a00aa1970,
    pop_city = a57ab1970,
    pop_urban = a57ac1970,
    pop_rural = a57ad1970,
    pop_white = b18aa1970,
    pop_black = b18ab1970,
    pop_aian = b18ac1970,
    pop_asian = b18ad1970,
    # pop_other = b18ae1970, # all NA
    pop_hisp = a35aa1970,
    vote_20 = state %in% c("Georgia", "Kentucky"), # can vote at 20
    # add to `cvap` to get better estimate
    cvap_adj = if_else(vote_20, a38au1970, -a38av1970 - a38aw1970),
    vap = rowSums(pick(a38ax1970:a38dw1970)) +
      if_else(vote_20, a38au1970 + a38av1970 + a38aw1970, 0),
    vap_white = rowSums(pick(aa1ax1970:aa1dw1970)) +
      if_else(vote_20, aa1au1970 + aa1av1970 + aa1aw1970, 0),
    vap_black = rowSums(pick(ac5ax1970:ac5dw1970)) +
      if_else(vote_20, ac5au1970 + ac5av1970 + ac5aw1970, 0)
  ) |>
  select(-vote_20)

c1970 <- c1970_raw3 |>
  left_join(c1970_raw1, by = c("statea", "countya")) |>
  left_join(c1970_raw2, by = c("statea", "countya")) |>
  filter(!statea %in% c("11", "02", "15")) |> # DC, AK, HI
  mutate(
    fips = str_c(statea, countya),
    .after = countya
  ) |>
  mutate( # fix Adams County, WI
    pop_city = coalesce(pop_city, 0),
    pop_urban = coalesce(pop_urban, 0),
    pop_rural = coalesce(pop_rural, pop),
  ) |>
  # Allocate this sliver to Park County
  rows_update(tibble(county = "Yellowstone National Park", fips = "30067"), by = "county") |>
  select(-statea, -countya) |>
  group_by(state, fips) |>
  summarize(
    county = county[1],
    across(-county, sum), .groups = "drop"
  )

# need to impute vap_black for some counties
m_imp <- lm(vap_black / vap ~ pop_black / pop, data = c1970)

c1970 <- c1970 |>
  mutate(
    vap_black = coalesce(vap_black, predict(m_imp, c1970)),
    cvap_unadj = cvap,
    cvap = cvap + cvap_adj,
    cvap_denom = cvap_denom + cvap_adj,
    cvap_black = coalesce(cvap_black, vap_black * cvap / cvap_denom),
    cvap_unadj = pmax(cvap_unadj, cvap_black + cvap_white),
    across(pop_city:pop_hisp, ~ . / pop),
    across(vap_white:vap_black, ~ . / vap),
    across(farm:nonfarm, ~ . / (farm + nonfarm)),
    across(educ_elem:educ_coll, ~ . / (
      educ_elem + educ_hsch + educ_coll)),
    across(cvap_white:cvap_black, ~ . / cvap_unadj),
    cvap = vap * (cvap / cvap_denom),
    across(inc_00_03k:inc_25_99k, ~ . / (
      inc_00_03k + inc_03_08k + inc_08_25k + inc_25_99k)),
  ) |>
  select(-cvap_adj, -cvap_denom, -cvap_unadj)

# sanity checks
with(c1970, {
  stopifnot(all(cvap <= vap))
  stopifnot(all(vap < pop))
  stopifnot(all(as.matrix(select(c1970, pop:inc_25_99k)) >= 0))
})

# https://digital.newberry.org/ahcb/downloads/gis/US_AtlasHCB_StateTerr_Gen05.zip
# https://digital.newberry.org/ahcb/downloads/gis/US_AtlasHCB_Counties_Gen05.zip
cty_shp <- read_sf("data-raw/US_AtlasHCB_Counties_Gen05/US_HistCounties_Gen05_Shapefile/") |>
  janitor::clean_names() |>
  filter(
    !state_terr %in% c("Hawaii", "District of Columbia", "Alaska"),
    !is.na(fips),
    # use Census Day; have checked that the 14 changes between 1968
    # election and then are not substantial
    end_date >= ymd("1970-04-01"),
    start_date <= ymd("1970-04-01")
  ) |>
  rmapshaper::ms_simplify(keep = 0.3, keep_shapes = TRUE) |>
  st_transform(5070) |>
  as_tibble() |>
  rows_update(tibble(id = "mts_yellowstonenp", fips = "30067"), by = "id") |>
  st_as_sf() |>
  group_by(fips) |>
  summarize(is_coverage = TRUE) |>
  as_tibble()

stopifnot(n_distinct(cty_shp$fips) == nrow(cty_shp))

c1970 <- full_join(c1970, cty_shp, by = "fips", relationship = "one-to-one") |>
  st_as_sf()

# https://www.icpsr.umich.edu/web/ICPSR/studies/8611
# for recode see first pages of codebook
icpsr_cty_recode <- tribble(
  ~state, ~county_fp, ~county_fp_override,
  "MARYLAND", "007", "009",
  "MARYLAND", "009", "011",
  "MARYLAND", "011", "013",
  "MARYLAND", "013", "015",
  "MARYLAND", "015", "017",
  "MARYLAND", "017", "019",
  "MARYLAND", "019", "021",
  "MARYLAND", "021", "023",
  "MARYLAND", "023", "025",
  "MARYLAND", "025", "027",
  "MARYLAND", "027", "029",
  "MARYLAND", "029", "031",
  "MARYLAND", "031", "033",
  "MARYLAND", "033", "035",
  "MARYLAND", "035", "039",
  "MARYLAND", "039", "041",
  "MARYLAND", "041", "043",
  "MARYLAND", "043", "045",
  "MARYLAND", "045", "047",
  "MISSOURI", "193", "186", # codebook has typo: 093 should be 193
  "NEVADA", "025", "510", # map Ormsby to Carson City
)
pres_cty <- haven::read_dta("data-raw/ICPSR_08611/DS0001/08611-0001-Data.dta",
  col_select = c(1:3, 664:669)
) |>
  transmute(
    state = haven::as_factor(V1) |>
      as.character() |>
      word(1, -2),
    county = V2,
    icpsr_code = V3,
    county_fp = if_else(str_sub(icpsr_code, -1) == "0",
      str_sub(icpsr_code, 1, -2), NA_character_
    ) |>
      as.integer() |>
      str_pad(3, "left", "0"),
    # factor of 10 error in MS
    pres_dem_hum = na_if(V664, 99.99) / if_else(state != "MISSISSIPPI", 100, 10),
    pres_rep_nix = na_if(V665, 99.99) / if_else(state != "MISSISSIPPI", 100, 10),
    pres_ind_wal = na_if(V666, 999.9) / 100,
    pres_oth = na_if(V667, 99.99) / 100,
    pres_total = na_if(V668, 9999999),
    pres_turn = na_if(V669, 999.9) / 100
  ) |>
  mutate(
    across(pres_dem_hum:pres_turn, ~ case_when(
      abs(. - 0.9999) < 1e-6 ~ NA_real_,
      abs(. - 9.999) < 1e-6 ~ NA_real_,
      TRUE ~ .
    )),
    norm_tot = coalesce(pres_dem_hum, 0) + coalesce(pres_rep_nix, 0) +
      coalesce(pres_ind_wal, 0) + coalesce(pres_oth, 0),
    across(pres_dem_hum:pres_oth, ~ . / norm_tot)
  ) |>
  select(-norm_tot) |>
  filter(
    icpsr_code != "9999", !is.na(pres_total),
    !state %in% c("ALASKA", "HAWAII")
  ) |>
  left_join(icpsr_cty_recode, by = c("state", "county_fp")) |>
  mutate(county_fp = if_else(is.na(county_fp_override), county_fp, county_fp_override)) |>
  select(-county_fp_override)

c1970 <- c1970 |>
  mutate(
    state_join = str_to_upper(state),
    county_fp = str_sub(fips, 3)
  )

st_info <- tibble::tribble(
  ~state, ~abbr, ~region, ~division,
  "Connecticut", "CT", "Northeast", "New England",
  "Maine", "ME", "Northeast", "New England",
  "Massachusetts", "MA", "Northeast", "New England",
  "New Hampshire", "NH", "Northeast", "New England",
  "Rhode Island", "RI", "Northeast", "New England",
  "Vermont", "VT", "Northeast", "New England",
  "Delaware", "DE", "South", "South Atlantic",
  "New Jersey", "NJ", "Northeast", "Middle Atlantic",
  "New York", "NY", "Northeast", "Middle Atlantic",
  "Pennsylvania", "PA", "Northeast", "Middle Atlantic",
  "Illinois", "IL", "North Central", "East North Central",
  "Indiana", "IN", "North Central", "East North Central",
  "Michigan", "MI", "North Central", "East North Central",
  "Ohio", "OH", "North Central", "East North Central",
  "Wisconsin", "WI", "North Central", "East North Central",
  "Iowa", "IA", "North Central", "West North Central",
  "Kansas", "KS", "North Central", "West North Central",
  "Minnesota", "MN", "North Central", "West North Central",
  "Missouri", "MO", "North Central", "West North Central",
  "Nebraska", "NE", "North Central", "West North Central",
  "North Dakota", "ND", "North Central", "West North Central",
  "South Dakota", "SD", "North Central", "West North Central",
  "Virginia", "VA", "South", "South Atlantic",
  "Alabama", "AL", "South", "East South Central",
  "Arkansas", "AR", "South", "West South Central",
  "Florida", "FL", "South", "South Atlantic",
  "Georgia", "GA", "South", "South Atlantic",
  "Louisiana", "LA", "South", "West South Central",
  "Mississippi", "MS", "South", "East South Central",
  "North Carolina", "NC", "South", "South Atlantic",
  "South Carolina", "SC", "South", "South Atlantic",
  "Texas", "TX", "South", "West South Central",
  "Kentucky", "KY", "South", "East South Central",
  "Maryland", "MD", "South", "South Atlantic",
  "Oklahoma", "OK", "South", "West South Central",
  "Tennessee", "TN", "South", "East South Central",
  "West Virginia", "WV", "South", "South Atlantic",
  "Arizona", "AZ", "West", "Mountain",
  "Colorado", "CO", "West", "Mountain",
  "Idaho", "ID", "West", "Mountain",
  "Montana", "MT", "West", "Mountain",
  "Nevada", "NV", "West", "Mountain",
  "New Mexico", "NM", "West", "Mountain",
  "Utah", "UT", "West", "Mountain",
  "Wyoming", "WY", "West", "Mountain",
  "California", "CA", "West", "Pacific",
  "Oregon", "OR", "West", "Pacific",
  "Washington", "WA", "West", "Pacific"
)

c1970 <- left_join(c1970, pres_cty,
  by = c("state_join" = "state", "county_fp"),
  relationship = "one-to-one", suffix = c("", "_pres")
) |>
  left_join(st_info, by = "state") |>
  select(-county_pres, -county_fp, -state_join, -icpsr_code) |>
  rename(pres_turn_icpsr = pres_turn) |>
  mutate(vap_other = 1 - vap_white - vap_black, .after = vap_black) |>
  mutate(cvap_other = 1 - cvap_white - cvap_black, .after = cvap_black) |>
  mutate(
    pres_abs = 1 - pres_dem_hum - pres_rep_nix - pres_ind_wal,
    .after = pres_ind_wal
  ) |>
  mutate(pres_turn = pres_total / cvap, .before = pres_turn_icpsr) |>
  relocate(geometry, .after = everything()) |>
  relocate(fips, .before = everything()) |>
  relocate(abbr, region, division, .after = state)

# check against official tally (after adding AK + HI + DC)
stopifnot(abs((sum(c1970$pres_total) + 83035 + 236218 + 170578) / 73199998 - 1) <= 1e-4)

write_rds(c1970, "data-raw/elec_1968.rds", compress = "xz")
