QR_WLS1<-function(X,k){
  M<-X
  m<-nrow(M)
  n<-ncol(M)
  initial<-svd(M)
  U<-initial$u[1:nrow(initial$u),1:k]
  V<-initial$v[1:nrow(initial$v),1:k]
  alpha<-4.685
  #alpha<-3
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
   # sigma<-1.4826*median(abs(r-median(abs(r))))
    sigma<-runmed(abs(RR),61)
    sigma[(N-30):N]<-calc.ends(abs(RR),30)
    plot(abs(RR), type='l')
    lines(sigma, type='l',col='red')
    #sigma<-6*median(abs(r))
    RR.trmatrix<-hankel(RR,L=120)
    #RR.row<-as.vector(RR.trmatrix)
    sigma.trmatrix<-hankel(sigma,L=120)
    #sigma.row<-as.vector(sigma.trmatrix)
   # plot(abs(RR.row),type='l')
   #lines(sigma.row,type='l',col='red')
   # plot(abs(RR.row)/sigma.row, type='l')
    plot(hankL1(abs(RR.trmatrix)/sigma.trmatrix), type='l')
    W<-weights1(RR.trmatrix,sigma.trmatrix,alpha,m,n)
    WW<-hankL1(W)
    plot(WW, type='l')
    Weights<-WW
    plot(Weights, type='l')
    title(main='Weights')
    
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


weights<-function(R,sigma,a,m,n){
  x<-R/sigma
  W<-matrix(0L, nrow = m, ncol = n)
  for (i in (1:nrow(W))){
    for (j in (1:ncol(W))){
      if (abs(x[i,j])<=a){
        W[i,j]<-(1-(abs(x[i,j])/a)^2)^2
      }
    }
  }
  #W[abs(x)<=a]<-(1-(abs(x)/a)^2)^2
  #W<-as.matrix(w,nrow=m,ncol=n,byrow=TRUE)
  return(W)
}


weights1<-function(R,sigma1,a,m,n){
  x<-R/sigma1
  W<-matrix(0L, nrow = m, ncol = n)
  for (i in (1:nrow(W))){
    for (j in (1:ncol(W))){
      if (sqrt(abs(x[i,j]))<=a){
        W[i,j]<-(1-(sqrt(abs(x[i,j]))/a)^2)^2
      }
    }
  }
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