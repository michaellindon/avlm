test_that("optimal_g minimizes t_radius", {
  for (n in c(10, 20, 50, 10000)) {
    k <- 5
    if (n <= k) next
    alpha <- 0.05
    nu <- n - k
    g_star <- optimal_g(n, k, alpha)
    f_min <- t_radius(g_star, n, alpha, nu)

    upper <- n * alpha^(2 / nu) / (1 - alpha^(2 / nu))
    grid <- exp(seq(log(upper) - 20, log(upper), length.out = 2000))
    grid_values <- sapply(grid, t_radius, n = n, alpha = alpha, nu = nu)
    expect_lte(f_min, min(grid_values) + 1e-8)
  }
})

test_that("optimal_g can go below 1 for small samples", {
  # With n = 10 and two coefficients the minimizer is about 0.68; a lower bound of 1 missed it.
  g_star <- optimal_g(10, 2, 0.05)
  expect_lt(g_star, 1)
})

test_that("optimal_g accepts the cluster degrees of freedom", {
  G <- 200
  g_star <- optimal_g(G, alpha = 0.05, nu = G - 1)
  grid <- exp(seq(log(g_star) - 2, log(g_star) + 2, length.out = 400))
  expect_lte(t_radius(g_star, G, 0.05, G - 1),
             min(sapply(grid, t_radius, n = G, alpha = 0.05, nu = G - 1)) + 1e-8)
})

test_that("optimal_g validates alpha", {
  expect_error(optimal_g(100, 2, 0), "alpha")
  expect_error(optimal_g(100, 2, 1), "alpha")
})
