elec_1968 = readRDS("data-raw/elec_1968.rds")
elec_1968$pres_oth[is.na(elec_1968$pres_oth)] = 0
elec_1968$adj = geomander::adjacency(st_make_valid(elec_1968)) |>
    lapply(\(x) elec_1968$fips[x + 1L])
elec_1968 = elec_1968 |>
    sf::st_drop_geometry() |>
    lapply(unname) |>
    tibble::as_tibble() |>
    dplyr::filter(abbr %in% c("LA", "AR", "TN", "MS", "AL", "GA", "FL", "SC", "NC", "VA", "TX"))

elec_1968 = elec_1968 |>
    ei_proportions(pop_white:pop_asian, .total=NULL) |>
    ei_proportions(vap_white:vap_other, .total=NULL) |>
    ei_proportions(cvap_white:cvap_other, .total=NULL)

usethis::use_data(elec_1968, overwrite = TRUE, compress="xz")

