devtools::load_all(".")
library(tidyverse)

n = 2000
p = 3
r2 = 0.25
specs = map(1:500, function(i) {
     ei_synthetic(
        n,
        p,
        n_x = n_x,
        r2_xz = r2,
        r2_bz = r2,
        b_cov = 0.0004 * (1 + diag(n_x))
    )
}, .progress = TRUE)

run_test = function(spec0, well_spec = TRUE, ...) {
    spec = ei_spec(
        spec0,
        predictors = starts_with("x"),
        outcome = starts_with("y"),
        total = attr(spec0, "ei_n"),
        covariates = if (well_spec) starts_with("z") else "z1",
        preproc = function(z) {
            if (ncol(z) == 0) {
                matrix(nrow = nrow(z), ncol = 0)
            } else {
                bases::b_tpsob(z, p = 25)
            }
        }
    )

    ei_test_car(spec, ...)
}

res = map(specs, ~ run_test(.x, well_spec = F, use_hc=F, iter=50), .progress = TRUE) |>
    bind_rows()
res2 = map(specs, ~ run_test(.x, well_spec = T, use_hc=F, iter=50), .progress = TRUE) |>
    bind_rows()

hist(res$p.value, breaks=50)
hist(res2$p.value, breaks=50)
mean(res$p.value <= 0.05)
mean(res2$p.value <= 0.05)
hist(pchisq(res$W, res$df, lower.tail=F), breaks=50)
plot(pchisq(res$W, res$df, lower.tail=F), res$p.value)
plot(res$p.value, res2$p.value)
plot(res$W, res2$W)
