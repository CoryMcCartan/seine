library(tidyverse)
library(rvest)
library(jsonlite)
library(here)

pg <- read_html("https://www.prri.org/research/2020-census-of-american-religion/")
urls <- pg |>
    html_elements("iframe[id*=datawrapper]") |>
    html_attr("src")
titles <- pg |>
    html_elements("iframe[id*=datawrapper]") |>
    html_attr("title") |>
    str_remove(", by County$") |>
    str_remove(" Identity$")
# remove overall
skip = str_detect(titles, "Religious Diversity Index")
titles = titles[!skip]
urls = urls[!skip]

pg <- read_html(urls[1])
js = pg |>
    html_element("script:not([type]):not([src])") |>
    html_text()
data_line = str_split_1(js, fixed("\n"))[1]
idx = str_locate(data_line, "parse\\(")[, "end"]
data_str = str_sub(data_line, idx + 1, -3)
d = fromJSON(fromJSON(data_str))$data$chartData |>
    read_delim(";", show_col_types=FALSE)

fn_renamer <- function(x) {
    x |>
        str_to_lower() |>
        str_replace("hispanic\\b", "hisp") |>
        str_replace("\\bevangelical protestant", "evang") |>
        str_replace("\\bprotestant", "prot") |>
        str_replace("\\bcatholic", "cath") |>
        str_replace("\\bchristian", "christ") |>
        str_remove(" mainline\\b") |>
        str_remove("^all ") |>
        str_replace("unaffiliated", "unaff") |>
        str_replace("latter-day saint \\(mormon\\)", "lds") |>
        str_replace_all(" ", "_") |>
        str_c("relig_", x=_)
}

d = d |>
    select(-Name, -Population) |>
    rename(fips=`FIPS-Code`) |>
    rename_with(fn_renamer, .cols=-fips) |>
    mutate(across(where(is.numeric), ~ . / 100)) |>
    select(-relig_white_christs)


write_csv(d, here("data-raw/cty_religion.csv"))
