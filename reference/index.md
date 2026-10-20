# Package index

## Anytime-valid models

Convert a fitted `lm` or `aov` model into an anytime-valid one.

- [`av()`](https://michaellindon.github.io/avlm/reference/av.md) :
  Anytime-valid Conversion Generic Function
- [`av(`*`<lm>`*`)`](https://michaellindon.github.io/avlm/reference/av.lm.md)
  : Convert a Linear Model (lm) Object to Anytime-Valid lm (avlm)
- [`av(`*`<aov>`*`)`](https://michaellindon.github.io/avlm/reference/av.aov.md)
  : Convert an aov Object to Anytime-Valid aov (avaov)

## Inference

- [`summary(`*`<avlm>`*`)`](https://michaellindon.github.io/avlm/reference/summary.avlm.md)
  : Summary Method for Anytime-Valid lm (avlm) Objects
- [`print(`*`<summary.avlm>`*`)`](https://michaellindon.github.io/avlm/reference/print.summary.avlm.md)
  : Print Method for summary.avlm Objects
- [`confint(`*`<avlm>`*`)`](https://michaellindon.github.io/avlm/reference/confint.avlm.md)
  : Confidence Intervals for Anytime-Valid lm (avlm) Objects
- [`vcov(`*`<avlm>`*`)`](https://michaellindon.github.io/avlm/reference/vcov.avlm.md)
  : Covariance Matrix for Anytime-Valid lm (avlm) Objects
- [`anova(`*`<avlm>`*`)`](https://michaellindon.github.io/avlm/reference/anova.avlm.md)
  : Anytime-Valid ANOVA Table for an avlm Object
- [`summary(`*`<avaov>`*`)`](https://michaellindon.github.io/avlm/reference/summary.avaov.md)
  : Summary Method for Anytime-Valid aov Objects
- [`print(`*`<summary.avaov>`*`)`](https://michaellindon.github.io/avlm/reference/print.summary.avaov.md)
  : Print Method for summary.avaov Objects

## Tuning

- [`optimal_g()`](https://michaellindon.github.io/avlm/reference/optimal_g.md)
  : Computes the value of g such that width of the \\1-\alpha\\
  confidence interval at sample size n is minimized
