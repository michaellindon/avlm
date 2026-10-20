## Submission

This is an update of avlm (0.1.0 -> 0.2.0). It adds cluster-robust anytime-valid inference
(`av(fit, cluster = ~ id)`), fixes the heteroskedasticity-robust F-test, and updates the
references in DESCRIPTION to the published JASA paper and the new arXiv paper on clustered
outcomes. See NEWS.md for details.

## R CMD check results

0 errors | 0 warnings | 0 notes

* local: macOS (aarch64-apple-darwin20), R 4.5.1, `R CMD check --as-cran`. The only NOTEs
  were local-environment ones (current time could not be verified; local HTML Tidy too old
  to validate the HTML manual), not issues with the package.

## Reverse dependencies

There are currently no reverse dependencies.
