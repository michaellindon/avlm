newton = Vectorize(function(alpha, n, p, phi, z2, x0, count=0, max_iter=10, tol=0.00000001){
  x1 = x0 - (log_bf(x0,n,p,phi,z2) + log(alpha) )/derivative_log_bf(x0,n,p,phi,z2)
  if(abs(x1-x0) < tol){
    return(x1)
  }else{
    if(count < max_iter){
      count = count + 1
      newton(alpha, n, p, phi, z2, x1, count = count)
    }else{
      return(x1)
    }
  }
})

projection_matrix = function(X){
  QRmat <- qr(X)
  Q <- qr.Q(QRmat)
  return(tcrossprod(Q))
}

MPxM = function(M, X){
  QRmat <- qr(X)
  Q <- qr.Q(QRmat)
  QtM = crossprod(Q,M)
  return(crossprod(QtM))
}

MP1M = function(M){
  n=dim(M)[1]
  Mbar = colSums(M)
  return(tcrossprod(Mbar)/n)
}