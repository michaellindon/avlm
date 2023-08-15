derivative_log_bf = function(t2, n, p, phi, z2){
  L = 1 / (n-p-1)
  U = L* phi / (phi + z2)
  return(0.5*(n-p)*(L/(1+t2*L) - U/(1+t2*U)))
}

log_bf = function(t2, n, p, phi, z2){
  0.5*log(phi / (phi + z2) ) + (0.5 * (n-p)) *(log(1 + t2 / (n-p-1)) - log(1 + phi * t2 / ((n-p-1) * (phi + z2 ))))
}

log_bf_multivariate = function(delta, n, p, d, Phi, ZtZ, s2){
  normalizing_constant =  if(d > 1) 0.5*log(det(Phi)) - 0.5*log(det(Phi+ZtZ)) else 0.5*log(detPhi) - 0.5*log(Phi+ZtZ)
  return (normalizing_constant +
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
