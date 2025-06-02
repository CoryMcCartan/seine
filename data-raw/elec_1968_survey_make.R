library(tidyverse)

# Data
datapath <- "data-raw/ICPSR_07508/DS0001/07508-0001-Data.txt"

# IDs -----
# Deck is the card number. case_id x state_id uniquely identifies person
card_id <- read_fwf(
    datapath,
    col_positions = fwf_cols(case_id = 3, state_id = 2,
                             region = 1, counting = 2, deck = 2)) |>
    mutate(row_id = 1:n(), .before = 1)


# Get to card level
card1s <- card_id |>
    filter(deck == "01") |>
    pull(row_id)

card3s <- card_id |>
    filter(deck == "03") |>
    pull(row_id)
card12s <- card_id |>
    filter(deck == "12") |>
    pull(row_id)

## Card specific variables
dat_card1 <- read_fwf(
    datapath,
    fwf_cols(case_id = 3,
             state_id = 2,
             intrastate_region = 1,
             counting = 2,
             deck = 2,
             sex = 1,
             age = 2,
             race = 1,
             nationality = 2,
             educ = 2,
             head_educ = 2,
             relig = 1
             )) |>
    slice(card1s)

dat_card3 <- read_fwf(
    datapath,
    fwf_cols(case_id = 3,
             state_id = 2,
             other3 = 47,
             voted68 = 1,
             novotereason68 = 1,
             pres68 = 1,
             intent68_novote = 1)) |>
    slice(card3s) |>
    select(-other3)

dat_card12 <- read_fwf(
    datapath,
    fwf_cols(case_id = 3,
             state_id = 2,
             other12 = 17,
             weight_intrastate = 4,
             weight_regional = 4,
             weight_national = 4
             )) |>
    slice(card12s) |>
    select(-other12) |>
    mutate(weight_national = as.numeric(weight_national),
           weight_intrastate = as.numeric(weight_intrastate))

# Join and recode -------
svy <- dat_card1 |>
    left_join(dat_card3, by = c("case_id", "state_id"),
              relationship = "one-to-one") |>
    left_join(dat_card12, by = c("case_id", "state_id"),
              relationship = "one-to-one") |>
    mutate(
        state = case_match(
            state_id,
            "02" ~ "Alabama",
            "04" ~ "Arizona",
            "05" ~ "Arkansas",
            "06" ~ "California",
            "07" ~ "Colorado",
            "08" ~ "Connecticut",
            "09" ~ "Delaware",
            "10" ~ "District of columbia",
            "11" ~ "Florida",
            "12" ~ "Georgia",
            "14" ~ "Idaho",
            "15" ~ "Illinois",
            "16" ~ "Indiana",
            "17" ~ "Iowa",
            "18" ~ "Kansas",
            "19" ~ "Kentucky",
            "20" ~ "Louisiana",
            "21" ~ "Maine",
            "22" ~ "Maryland",
            "23" ~ "Massachusetts",
            "24" ~ "Michigan",
            "25" ~ "Minnesota",
            "26" ~ "Mississippi",
            "27" ~ "Missouri",
            "28"  ~ "Montana",
            "29"  ~ "Nebraska",
            "30"  ~ "Nevada",
            "31"  ~ "New Hampshire",
            "32"  ~ "New Jersey",
            "33"  ~ "New Mexico",
            "34"  ~ "New York",
            "35"  ~ "North Carolina",
            "37"  ~ "Ohio",
            "38"  ~ "Oklahoma",
            "39"  ~ "Oregon",
            "40"  ~ "Pennsylvania",
            "41"  ~ "Rhode Island",
            "42"  ~ "South Carolina",
            "43"  ~ "South Dakota",
            "44"  ~ "Tennessee",
            "45"  ~ "Texas",
            "46"  ~ "Utah",
            "47"  ~ "Vermont",
            "48"  ~ "Virginia",
            "49"  ~ "Washington",
            "50"  ~ "West Virginia",
            "51"  ~ "Wisconsin",
            "52"  ~ "Wyoming"
        ),
        cfdr_states = state %in% c(
            "Alabama", "Mississippi", "North Carolina",
            "Georgia", "South Carolina", "Louisiana",
            "Texas", "Florida", "Arkansas", "Virginia",
            "Tennessee"
        ),
        sex = case_match(
            sex,
            1 ~ "Male",
            2 ~ "Female"
        ),
        educ = factor(
            educ,
           labels = c("First grade",
            "Second grade",
            "Third grade",
            "Fourth grade",
            "Fifth grade",
            "Sixth grade",
            "Seventh grade",
            "Eighth grade",
            "Ninth grade",
            "Tenth grade",
            "Eleventh grade",
            "Twelveth grade",
            "Freshman ",
            "Sophomore",
            "Junior",
            "Senior",
            "One year beyond college",
            "Two or more years beyond college",
            "No schooling",
            "Not sure how much schooling",
            "Not ascertained (Error)")
        ),
        race = case_match(
            race,
            1 ~ "White",
            2 ~ "Negro",
            3 ~ "Mexican-American",
            4 ~ "Puerto-Rican",
            5 ~ "Oriental",
            6 ~ "American Indian",
            7 ~ "Other (specify)",
            8 ~ "Not sure"
        ),
        voted68 = case_match(
            voted68,
            1 ~ "Voted",
            2 ~ "Did not vote",
            3 ~ "Refused"
        ),
        pres68_full = case_match(
            pres68,
            1 ~ "Nixon",
            2 ~ "Humphrey",
            3 ~ "Wallace",
            4 ~ "Other (vol)",
            5 ~ "Don't recall",
            6 ~ "Didn't vote for President",
            7 ~ "Refused",
            9 ~ "Not ascertained (Error)",
            0 ~ NA_character_),
        novotereason68 = case_match(
            novotereason68,
            1 ~ "Not registered;  not  eligible;  couldn't  meet residence requirements.",
            2 ~ "Illness (self or family); old age.",
            3 ~ "Too busy; didn't have time; out of town.",
            4 ~ "Never vote; general lethargy.",
            5 ~ "Didn't like candidates; no  difference  between candidates; didn't care who won.",
            6 ~ "Other",
            7 ~ "Not sure",
            9 ~ "Not ascertained (Error)",
            0 ~ NA_character_
        )
        ) |>
    # factors
    mutate(
        ## lump together the smalle candidates
        pres68 = fct_infreq(fct_lump_n(pres68_full, 4)),
        race = fct_infreq(race),
        voted68 = fct_infreq(voted68),
        age = as.numeric(age)) |>
    relocate(case_id, state_id, state)


# Tabulations -----
# unweighted counts
xtabs( ~ race + pres68, svy, subset = cfdr_states, addNA = TRUE)

# weighted proportions
xtabs(weight_national ~ race + pres68, svy, subset = cfdr_states) |>
    prop.table(1) |> round(3)

# unweighted
xtabs(~ race + voted68, svy, addNA = TRUE,
      subset = cfdr_states) |>
    addmargins()


# Export ----
svy |>
    select(state_id, case_id,
           state,
           cfdr_states,
           age, sex, educ,
           race,
           matches("68"),
           matches("weight")) |>
    write_rds("data-raw/elec_1968_svy.rds")

