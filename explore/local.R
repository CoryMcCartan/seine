devtools::load_all(".")
library(tidyverse)

data(elec_1968)
elec_1968$vap_nonwhite = 1 - elec_1968$vap_white

spec = ei_spec(elec_1968, c(vap_white, vap_nonwhite), pres_dem_hum,
               total = pres_total, covariates = c(state, pop_urban, farm))

m = ei_ridge(spec)

ei_est_local(m, spec, conf_level = 0.95)

idx = which.max(rl$preds$vap_white < 0)
y = spec$pres_dem_hum[idx]
x = rl$x[idx, ]
b0 = sapply(rl$preds, \(x) x[idx, ])
all.equal(rl$yhat[idx], sum(x * b))

fn_dist <- function(b, x, y) {
    eps = y - sum(x * b)
    1e3*(eps^2 + sum(b[b < 0]^2) + sum((b[b > 1] - 1)^2)) + sum((b - b0)^2)
}
gr_dist <- function(b, x, y) {
    eps = y - sum(x * b)
    1e3*(-2*eps*x + 2*b*(b < 0) + 2*(b-1)*(b > 1)) + 2*(b - b0)
}
(c(fn_dist(b0 + c(1, 0)*1e-6, x, y), fn_dist(b0 + c(0, 1)*1e-6, x, y)) - fn_dist(b0, x, y)) / 1e-6
gr_dist(b0, x, y)

b = optim(b0, fn_dist, gr_dist, x=x, y=y, method="L-BFGS-B")$par |> print()
fn_dist(b, x, y)

crossing(
    b0 = seq(-0.5, 1.5, 0.01),
    b1 = seq(-0.5, 1.5, 0.01),
) |>
    mutate(
        eps = y - b0*x[1] - b1*x[2],
        pen = eps^2 + if_else(b0 < 0, b0^2, 0) + if_else(b1 < 0, b1^2, 0) +
            if_else(b0 > 1, (b0-1)^2, 0) + if_else(b1 > 1, (b1-1)^2, 0)
    ) |>
ggplot(aes(b0, b1, z=sqrt(pen))) +
    geom_contour()
    geom_raster()



