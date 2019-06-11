QR_WLS<-function(X,k){
  M<-X
  m<-nrow(M)
  n<-ncol(M)
  initial<-svd(M)
  U<-initial$u[1:nrow(initial$u),1:k]
  V<-initial$v[1:nrow(initial$v),1:k]
  alpha<-4.685
  eps<-1e-8
  eps1<-1e-5
  maxITER<-20
  maxiter<-10
  ITER<-0
  iter<-0
  
  repeat {
    R<-M-U%*%t(V)
    r<-as.vector(t(R))
    sigma<-1.4826*median(abs(r-median(abs(r))))
    W<-weights(R,sigma,alpha,m,n)
    
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
  }
  M_est<-U%*%t(V)
  return(M_est)
}

# weights1<-function(R,sigma1,a,m,n){
#   x<-R/sigma1
#   W<-matrix(0L, nrow = m, ncol = n)
#   W[abs(x)<=a]<-(1-(abs(x)/a)^2)^2
#   return(W)
# }

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