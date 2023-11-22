log_bf = function(t2, nu, phi, z2){
  r <- phi / (phi + z2)
  0.5*log(r) + (0.5 * (nu+1)) * (log(1 + t2 / nu) - log(1 + r * t2 / nu  ))
}

log_bf_multivariate = function(delta, n, p, d, Phi, ZtZ, s2){
  return (0.5*log(det(Phi)) - 0.5*log(det(Phi+ZtZ)) + 
        (0.5*(n-p))*(
          log(1+t(delta) %*% ZtZ %*% delta / (s2 * (n-p-d))) - 
          log(1+t(delta) %*% (ZtZ - ZtZ %*% solve(ZtZ + Phi, ZtZ)) %*% delta / (s2 * (n-p-d)))))
}

multivariate_F_sequential_p_value = function(lmfit, phi, s2) {
  W = model.matrix(lmfit)
  p = 1
  n = dim(W)[1]
  d = dim(W)[2] - 1
  Z = W[, 2:(d + p)]
  Y = lmfit$residuals + lmfit$fitted.values
  delta = lmfit$coefficients[2:(d + p)]
  ZtZ = MPxM(Z,W) - MP1M(Z)
  Phi = diag(d) * phi
  return (min(1, exp(
    -1 * log_bf_multivariate(delta, n, p, d, Phi, ZtZ, s2)
  )))
}