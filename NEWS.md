# avlm 0.1.0

## New Features

- First release of the **avlm** package, which implements anytime-valid inference for linear models and ANOVA.
- Provides a drop-in replacement for `summary.lm()` and `summary.aov()` with sequential p-values and confidence intervals.
- Includes support for:
  - `av()`: applies the anytime-valid transformation to `lm` and `aov` objects.
  - `confint()`: computes anytime-valid confidence intervals for regression coefficients.

## Infrastructure

- Full documentation with `roxygen2`.
- Unit tests using `testthat`.

## Notes

- Based on the methodology described in  
  *"Anytime-Valid Linear Models and Regression Adjusted Causal Inference in Randomized Experiments"*  
  available at: <https://arxiv.org/abs/2210.08589>


