weights<-function(R,sigma,a,m,n){
  x<-R/sigma
  W<-matrix(0L, nrow = m, ncol = n)
  W[which(abs(x)<=a,arr.ind=TRUE)]<-(1-(abs(x[which(abs(x)<=a,arr.ind=TRUE)])/a)^2)^2
  return(W)
}

calc.ends<-function(vec,L){
  N<-length(vec)
  res<-rep(NaN, (L+1))
  for (i in (1:(L+1))) {
    slide<-vec[(N-2*L+i):N]
    res[i]<-median(slide)
  }
  return(res)
}


QR_WLS<-function(X,k,trend.ver='loess'){
  M<-X
  m<-nrow(M)
  n<-ncol(M)
  initial<-svd(M)
  U<-initial$u[1:nrow(initial$u),1:k]
  Lambda<-initial$d[1:k]
  V<-initial$v[1:nrow(initial$v),1:k]
  U<-U%*%diag(Lambda)
  alpha<-4.046
  eps<-1e-8
  eps1<-1e-5
  maxITER<-20
  maxiter<-10
  ITER<-0
  iter<-0
  
  
  repeat {
    R<-M-U%*%t(V)
    r<-as.vector(t(R))
    RR<-hankL1(R)
    loessMod30 <- loess(abs(RR) ~ c(1:length(abs(RR))), span=0.30)
    sigma <- predict(loessMod30)
    
    sigma1<-runmed(abs(RR),61)
    sigma1[(N-30):N]<-calc.ends(abs(RR),30)
    
    sigma2<-real.trend*sqrt(2/pi)
    
    RR.trmatrix<-hankel(RR,L=120)
    if (trend.ver == 'loess') {sigma.trmatrix<-hankel(sigma,L=120)}
    else if (trend.ver == 'median') {sigma.trmatrix<-hankel(sigma1,L=120)}
    else if (trend.ver == 'real') {sigma.trmatrix<-hankel(sigma2,L=120)}
    #RR.row<-as.vector(RR.trmatrix)
    #sigma.trmatrix<-hankel(sigma,L=120)
    #plot(hankL1(abs(RR.trmatrix)/sigma.trmatrix), type='l')
    W<-weights(RR.trmatrix,sigma.trmatrix,alpha,m,n)
    WW<-hankL1(W)
    Weights<-WW
    #plot(Weights, type='l')
    #title(main='Weights')
    
    repeat{
      for (i in (1:m)){
        Wi<-diag(W[i,1:ncol(W)])
        mi<-M[i,1:ncol(M)]
        QR <- qr(t(V)%*%Wi%*%V)
        Q <- qr.Q( QR )
        R <- qr.R( QR )
        beta<-backsolve(R,t(Q)%*%t(V)%*%Wi%*%mi)
        U[i,1:ncol(U)]<-beta
      }
      U<-U[1:nrow(U),1:k]
      
      for (j in (1:n)){
        Wj<-diag(W[1:nrow(W),j])
        mj<-M[1:nrow(M),j]
        QR <- qr(t(U)%*%Wj%*%U)
        Q <- qr.Q( QR )
        R <- qr.R( QR )
        beta<-backsolve(R,t(Q)%*%t(U)%*%Wj%*%mj)
        V[j,1:ncol(V)]<-beta
      }
      V<-V[1:nrow(V),1:k]
      
      iter<-iter+1
      if ( ((frobenius.norm(W^{1/2}*(M-U%*%t(V))))^2<eps) | (iter > maxiter) ) {break}
    }
    
    ITER<-ITER+1
    if ( ((frobenius.norm(W^{1/2}*(M-U%*%t(V))))^2<eps) | (ITER > maxITER) ) { break}
    iter<-0
  }
  M_est<-U%*%t(V)
  return(M_est)
}