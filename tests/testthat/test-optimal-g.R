library(testthat)

test_that("optimal_g returns a value that minimizes t_radius", {
  # Set up parameters: n, number_of_coefficients, and alpha.
  n <- 10000
  number_of_coefficients <- 5
  alpha <- 0.05
  nu <- n - number_of_coefficients
  
  # Compute the optimal g using the function
  g_star <- optimal_g(n, number_of_coefficients, alpha)
  
  # Evaluate t_radius at the optimal g
  f_min <- t_radius(g_star, nu, n, alpha)
  
  # Create a small grid around g_star.
  # We choose a relative shift (here 1% of g_star) to test local minimality.
  delta <- max(g_star * 0.01, 1e-4)
  g_left <- ifelse(g_star - delta > 1, g_star - delta, 1)  # ensure we don't go below 1 (our lower bound)
  g_right <- g_star + delta
  
  f_left <- t_radius(g_left, n, number_of_coefficients, alpha)
  f_right <- t_radius(g_right, n, number_of_coefficients, alpha)
  
  # Check that the value at g_star is less than or equal to values at adjacent grid points.
  expect_true(f_min <= f_left,
              info = sprintf("t_radius at g_star (%g) = %g should be <= t_radius at g_left (%g) = %g", 
                             g_star, f_min, g_left, f_left))
  expect_true(f_min <= f_right,
              info = sprintf("t_radius at g_star (%g) = %g should be <= t_radius at g_right (%g) = %g", 
                             g_star, f_min, g_right, f_right))
  
  # Optionally, check over a larger grid around g_star.
  grid_points <- seq(max(1, g_star - 5*delta), g_star + 5*delta, length.out = 25)
  grid_values <- sapply(grid_points, function(g) t_radius(g, n, number_of_coefficients, alpha))
  
  expect_true(all(grid_values >= f_min),
              info = "g_star does not provide a local minimum over a reasonable grid around it.")
})
