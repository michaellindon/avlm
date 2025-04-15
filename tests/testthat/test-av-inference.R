test_that("anytime-valid p-values are more conservative than standard p-values for linear models", {
  # Test with different precision parameters
  for (g_val in c(1, 2, 5)) {
    # Fit standard and anytime-valid models
    std_fit <- lm(mpg ~ wt + hp, data = mtcars)
    std_summary <- summary(std_fit)
    std_pvals <- std_summary$coefficients[, 4]
    
    av_fit <- av(std_fit, g = g_val)
    av_summary <- summary(av_fit)
    av_pvals <- av_summary$coefficients[, 4]
    
    # Check if all anytime-valid p-values are greater than or equal to standard p-values
    for (i in seq_along(std_pvals)) {
      expect_gte(av_pvals[i], std_pvals[i])
    }
  }
})

test_that("anytime-valid p-values are more conservative than standard p-values for ANOVA models", {
  # Test with different precision parameters
  for (g_val in c(1, 2, 5)) {
    # Fit standard and anytime-valid models
    std_fit <- aov(Sepal.Length ~ Species, data = iris)
    std_summary <- summary(std_fit)
    std_pvals <- std_summary[[1]]$`Pr(>F)`[1]  # Extract p-value for Species
    
    av_fit <- av(std_fit, g = g_val)
    av_summary <- summary(av_fit)
    av_pvals <- av_summary[[1]]$`Pr(>F)`[1]  # Extract p-value for Species
    
    # Check if anytime-valid p-value is greater than or equal to standard p-value
    expect_gte(av_pvals, std_pvals)
  }
})

test_that("anytime-valid confidence intervals are supersets of standard confidence intervals for linear models", {
  # Test with different precision parameters and confidence levels
  for (g_val in c(1, 2, 5)) {
    for (level_val in c(0.90, 0.95, 0.99)) {
      # Fit standard and anytime-valid models
      std_fit <- lm(mpg ~ wt + hp, data = mtcars)
      std_ci <- confint(std_fit, level = level_val)
      
      av_fit <- av(std_fit, g = g_val)
      av_ci <- confint(av_fit, level = level_val)
      
      # Check if all anytime-valid CIs contain standard CIs
      for (i in 1:nrow(std_ci)) {
        coef_name <- rownames(std_ci)[i]
        
        # Lower bound should be smaller or equal
        expect_lte(av_ci[coef_name, 1], std_ci[coef_name, 1])
        
        # Upper bound should be greater or equal
        expect_gte(av_ci[coef_name, 2], std_ci[coef_name, 2])
      }
    }
  }
})

test_that("anytime-valid model classes are correctly assigned", {
  # Linear model
  std_lm <- lm(mpg ~ wt + hp, data = mtcars)
  av_lm <- av(std_lm)
  expect_s3_class(av_lm, "avlm")
  expect_true("lm" %in% class(av_lm))
  
  # ANOVA model
  std_aov <- aov(Sepal.Length ~ Species, data = iris)
  av_aov <- av(std_aov)
  expect_s3_class(av_aov, "avaov")
  expect_true("aov" %in% class(av_aov))
})

test_that("g parameter is correctly stored and accessed", {
  # Test different g values
  for (g_val in c(1, 2, 5)) {
    # Linear model
    std_lm <- lm(mpg ~ wt + hp, data = mtcars)
    av_lm <- av(std_lm, g = g_val)
    expect_equal(attr(av_lm, "g"), g_val)
    
    # ANOVA model
    std_aov <- aov(Sepal.Length ~ Species, data = iris)
    av_aov <- av(std_aov, g = g_val)
    expect_equal(attr(av_aov, "g"), g_val)
  }
})

test_that("summary method retains anytime-valid properties", {
  # Linear model
  std_lm <- lm(mpg ~ wt + hp, data = mtcars)
  av_lm <- av(std_lm, g = 2)
  av_summary <- summary(av_lm)
  expect_s3_class(av_summary, "summary.avlm")
  expect_equal(attr(av_summary, "g"), 2)
  
  # ANOVA model
  std_aov <- aov(Sepal.Length ~ Species, data = iris)
  av_aov <- av(std_aov, g = 3)
  av_summary <- summary(av_aov)
  expect_s3_class(av_summary, "summary.avaov")
  expect_equal(attr(av_summary, "g"), 3)
})

test_that("anytime-valid methods work with increasingly complex models", {
  # More complex linear model with interactions
  complex_lm <- lm(mpg ~ wt * hp + factor(cyl), data = mtcars)
  av_complex_lm <- av(complex_lm)
  expect_no_error(summary(av_complex_lm))
  expect_no_error(confint(av_complex_lm))
  
  # More complex ANOVA model
  complex_aov <- aov(Sepal.Length ~ Species * Petal.Width, data = iris)
  av_complex_aov <- av(complex_aov)
  expect_no_error(summary(av_complex_aov))
})