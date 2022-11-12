test_that("standard linear algebra is equal to faster linear algebra", {
  set.seed(1)
  n=100
  X = replicate(3, rnorm(n))
  trt = sample(c(0,1), size = n, replace = TRUE)
  y = 1 + 2*X[,1] + 3*X[, 2]+ 4*X[, 3]+ 2.3*trt + 10*rnorm(n)
  df = data.frame(X)
  df$trt = trt
  df$y = y
  lmfit = lm(y~.,data=df)
  
  W = model.matrix(lmfit)
  p = 1
  n = dim(W)[1]
  d = dim(W)[2] - 1
  Z = W[, 2:(d + p)]
  Y = lmfit$residuals + lmfit$fitted.values
  delta = lmfit$coefficients[2:(d + p)]
  
  Pw = W %*% solve(t(W) %*% W, t(W))
  Px = matrix(1, n, n) / n
  ZtZ1 = t(Z) %*% (Pw - Px) %*% Z
  
  ZtZ2 = MPxM(Z,W) - MP1M(Z)
  
  expect_equal(ZtZ1, ZtZ2)
})

test_that("Projection Matrix", {
  set.seed(1)
  n=100
  W = replicate(3, rnorm(n))
  Pw1 = W %*% solve(t(W) %*% W, t(W))
  Pw2 = projection_matrix(W)
  expect_equal(Pw1,Pw2)
})

test_that("Symmetric Projection Matrix Product",{
  set.seed(1)
  n=100
  W = replicate(3, rnorm(n))
  Z = replicate(3, rnorm(n))
  Pw = projection_matrix(W)
  out1 = t(Z)%*%Pw%*%Z
  out2 = MPxM(Z,W)
  expect_equal(out1, out2)
})

test_that("Fast", {
  set.seed(1)
  n=100
  W = replicate(3, rnorm(n))
  P1 = matrix(1, n, n) / n
  
  out1 = MP1M(W)
  out2 = t(W) %*% P1 %*% W
  expect_equal(out1,out2)
})