## CRAN submission of avlm 0.1.0

This is the first submission of the **avlm** package to CRAN.

The package implements anytime-valid inference for linear models and ANOVA, as described in:
  Lindon (2022), "Anytime-Valid Linear Models and Regression Adjusted Causal Inference in Randomized Experiments"
  https://arxiv.org/abs/2210.08589

### Local checks

- [x] R CMD check passed with `--as-cran` on macOS (no ERRORs, WARNINGs, or NOTEs)
- [x] devtools::check_win_devel() passes
- [x] rhub::check_for_cran() passes on Linux and Windows

### Comments

- This is a new package with no reverse dependencies.
- No external software or system dependencies are required.

