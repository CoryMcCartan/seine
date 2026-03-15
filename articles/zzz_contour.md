# Custom sensitivity contour plots

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
objects returned by
[`ei_sens()`](https://corymccartan.com/seine/reference/ei_sens.md)
produces a standard sensitivity contour plot. For publication-quality
figures with more control over styling, labeled contour lines, and
benchmark annotations, a custom `ggplot2` implementation is often
preferable. This vignette shows a generic version of such a plot, using
the [`geomtextpath`](https://cran.r-project.org/package=geomtextpath)
package for labeled contours and
[`ggrepel`](https://cran.r-project.org/package=ggrepel) for
non-overlapping benchmark labels.

The code below assumes `sens` is an object returned by
[`ei_sens()`](https://corymccartan.com/seine/reference/ei_sens.md) and
`bench` is a benchmarking data frame from
[`ei_bench()`](https://corymccartan.com/seine/reference/ei_bench.md),
both filtered to a single outcome of interest. Both are described in
[`vignette("sensitivity")`](https://corymccartan.com/seine/articles/sensitivity.md).

``` r
library(ggplot2)
library(geomtextpath)
library(ggrepel)

# --- filter to one outcome --------------------------------------------------
outcome_name = "my_outcome"
d_sens = subset(sens, outcome == outcome_name)
d_bench = subset(bench, outcome == outcome_name)

# clamp any Inf bias values for plotting
d_sens$bias_bound = with(d_sens, ifelse(
    is.finite(bias_bound),
    bias_bound,
    1.5 * max(bias_bound[is.finite(bias_bound)])
))

# --- contour break structure -------------------------------------------------
# major breaks get text labels; semi-major and minor add visual density
contour_exp = -1:1
breaks_minor   = 10^c(contour_exp, 2) %x% c(2:4, 6:9)
breaks_semimaj = 10^contour_exp %x% 5
breaks_maj     = 10^contour_exp %x% 10
c_col = "#f08f88"   # contour line color

# --- build plot --------------------------------------------------------------
ggplot(d_sens, aes(c_predictor, c_outcome)) +
    # minor contour lines (no labels)
    geom_contour(aes(z = bias_bound), breaks = breaks_minor,
                 color = c_col, lwd = 0.12) +
    # semi-major contour lines (no labels)
    geom_contour(aes(z = bias_bound), breaks = breaks_semimaj,
                 color = c_col, lwd = 0.48) +
    # major contour lines with numeric labels
    geom_textcontour(
        aes(z = bias_bound),
        breaks = breaks_maj,
        color = c_col, lwd = 0.7,
        hjust = 0.55, vjust = 1.25,
        halign = "left", size = 3.0,
        fontface = "bold", upright = TRUE, straight = TRUE
    ) +
    # a labeled contour for a specific reference value (e.g. the point estimate)
    geom_textcontour(
        aes(z = bias_bound, label = "Estimated effect"),
        breaks = d_sens$estimate[1],
        color = "#425682", lwd = 0.7,
        hjust = 0.55, vjust = 1.25,
        size = 3.0, fontface = "bold"
    ) +
    # benchmark points and labels
    geom_point(data = d_bench, size = 0.7) +
    geom_text_repel(
        aes(label = covariate),
        data = d_bench,
        hjust = 1.25, nudge_x = 0.002, nudge_y = 0.002,
        size = 3
    ) +
    # square-root scale spreads out the lower-left corner
    scale_x_continuous(breaks = seq(0, 1, 0.1), transform = "sqrt") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), transform = "sqrt") +
    coord_cartesian(xlim = 0:1, ylim = 0:1, expand = FALSE) +
    labs(
        x = bquote(1 - {R^2}[alpha^A ~ "~" ~ alpha]),
        y = bquote({R^2}[bar(Y) ~ "~" ~ A ~ "|" ~ bar(X) ~ "," ~ Z]),
        title = paste("Sensitivity contour plot:", outcome_name)
    ) +
    theme_bw() +
    theme(
        panel.grid = element_line(color = "#aaa", linetype = "dotted",
                                  linewidth = 0.24),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12)
    )
```
