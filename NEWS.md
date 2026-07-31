# seine 0.1.2.9999

* Add `ei_target()`, which adjusts a fitted regression using Targeted
  Minimum-Loss Estimation (TMLE), an alternative to the DML one-step estimator 
  that incorporates bounds and constraints on the outcome.
  The resulting regression function can be provided to `ei_est()` to produce
  global estimates that respect bounds, and can be provided to `ei_est_local()`
  to incorporate some double robustness.
* Support targeted models as well as external regression models passed through
  `ei_wrap_model()` in the sensitivity analyses
* Fix bug in sensitivity benchmarking parameter calculation

# seine 0.1.2

* Double/debiased ML for ecological inference
* Sensitivity analysis for EI
* Partial identification bounds
* Local (unit-level) estimates and confidence intervals
* Preprocessing helper functions
* Synthetic EI data
