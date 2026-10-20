test_that("avaov objects give anytime-valid confidence intervals and covariances", {
  fit <- aov(Sepal.Length ~ Species, data = iris)
  av_fit <- av(fit, g = 2)
  expect_equal(confint(av_fit), confint(av(lm(Sepal.Length ~ Species, data = iris), g = 2)))
  expect_false(isTRUE(all.equal(confint(av_fit), confint(fit))))
  expect_equal(vcov(av_fit), vcov(fit))
})

test_that("anova() on an avlm object replaces the p-values with anytime-valid ones", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  av_fit <- av(fit, g = 3)
  tab <- anova(av_fit)
  classical <- anova(fit)
  expect_equal(tab$`F value`, classical$`F value`)
  expect_equal(tab$Df, classical$Df)
  # same e-values as summary.avaov
  expected <- summary(av(aov(mpg ~ wt + hp, data = mtcars), g = 3))[[1]]$`Pr(>F)`
  expect_equal(tab$`Pr(>F)`, expected)
  expect_true(all(tab$`Pr(>F)`[1:2] >= classical$`Pr(>F)`[1:2]))
  expect_true(is.na(tab$`Pr(>F)`[3]))
  expect_output(print(tab), "Anytime-valid")
  expect_error(anova(av_fit, fit), "single model")
})

test_that("avaov summaries report g and reject multivariate responses", {
  expect_output(print(summary(av(aov(Sepal.Length ~ Species, data = iris), g = 2.5))), "g = 2.5")
  expect_error(av(aov(cbind(mpg, qsec) ~ factor(cyl), data = mtcars)), "multivariate")
})

test_that("too few clusters for the coefficients gives a warning", {
  fit <- lm(mpg ~ wt + hp + qsec, data = mtcars)
  expect_warning(summary(av(fit, cluster = ~ am)), "singular")
  expect_warning(confint(av(fit, cluster = ~ am)), "singular")
  expect_no_warning(summary(av(lm(mpg ~ wt, data = mtcars), cluster = ~ cyl)))
})

test_that("confint validates level", {
  expect_error(confint(av(lm(mpg ~ wt, data = mtcars)), level = 1.5), "level")
})
