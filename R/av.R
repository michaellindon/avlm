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
  M = if(is.null(n)) matrix(M, ncol=1) else M
  Mbar = colSums(M)
  return(tcrossprod(Mbar)/n)
}
