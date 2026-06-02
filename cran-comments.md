## Test environments
* local R installation (macOS), R 4.6.0
* win-builder (devel)
* windows-latest (on gh-actions), (release)
* macos-latest (on gh-actions), (release)
* ubuntu-latest (on gh-actions), (release)
* ubuntu-latest (on gh-actions), (devel)
* ubuntu-latest (on gh-actions), (old-release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a resubmission of a new package which moves `par()` inside `plot.ei_sens()`
  from the end of the function to an `on.exit()` call right before plotting 
  parameters are changed.
