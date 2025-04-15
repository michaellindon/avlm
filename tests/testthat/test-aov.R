test_that("anytime-valid p-values are more conservative than standard p-values for ANOVA models", {
  # Test with different precision parameters
  for (g_val in c(1, 2, 5)) {
    # Fit standard and anytime-valid ANOVA models
    std_fit <- aov(Sepal.Length ~ Species, data = iris)
    std_summary <- summary(std_fit)
    std_pvals <- std_summary[[1]]$`Pr(>F)`[1]  # Extract the p-value for Species
    
    av_fit <- av(std_fit, g = g_val)
    av_summary <- summary(av_fit)
    av_pvals <- av_summary[[1]]$`Pr(>F)`[1]
    
    # Check if the anytime-valid p-value is greater than or equal to the standard p-value
    expect_gte(av_pvals, std_pvals)
  }
})

test_that("anytime-valid methods work with increasingly complex ANOVA models", {
  # Test using a more complex ANOVA model
  complex_aov <- aov(Sepal.Length ~ Species * Petal.Width, data = iris)
  av_complex_aov <- av(complex_aov)
  expect_no_error(summary(av_complex_aov))
})
