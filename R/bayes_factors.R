log_G_t = function(t2, nu, n, g){
  r <- g / (g + n)
  0.5*log(r) + (0.5 * (nu+1)) * (log(1 + t2 / nu) - log(1 + r * t2 / nu  ))
}

p_G_t = function(log_G_t_values){
  return(min(1.0, exp(-log_G_t_values)))
}

log_G_f = function(f, d, nu, n, g){
  r <- g / (g + n)
  0.5 * d * log(r) + (0.5 * (nu+d)) * (log(1 + (d / nu) * f) - log(1 + r * (d / nu) * f  ))
}

p_G_f = function(log_G_f_values){
  return(min(1.0, exp(-log_G_f_values)))
}

p_t = function(t2, nu, phi, z2){
  pmin(1, exp(-1 * log_E_t(t2, nu, phi, z2)))
}


log_E_t = function(t2, nu, phi, z2){
  r <- phi / (phi + z2)
  0.5*log(r) + (0.5 * (nu+1)) * (log(1 + t2 / nu) - log(1 + r * t2 / nu  ))
}

p_t = function(t2, nu, phi, z2){
  pmin(1, exp(-1 * log_E_t(t2, nu, phi, z2)))
}

log_E_f = function(delta, n, p, d, Phi, ZtZ, s2){
  normalizing_constant =  if(d > 1) 0.5*log(det(Phi)) - 0.5*log(det(Phi+ZtZ)) else 0.5*log(Phi) - 0.5*log(Phi+ZtZ)
  sol = if(d>1) solve(ZtZ + Phi, ZtZ) else ZtZ / (ZtZ+Phi)
  return (normalizing_constant +
        (0.5*(n-p))*(
          log(1+t(delta) %*% ZtZ %*% delta / (s2 * (n-p-d))) - 
          log(1+t(delta) %*% (ZtZ - ZtZ %*% sol) %*% delta / (s2 * (n-p-d)))))
}

p_F = function(lmfit, phi, s2) {
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
    -1 * log_E_f(delta, n, p, d, Phi, ZtZ, s2)
  )))
}
