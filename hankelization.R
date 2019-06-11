#L1-hankelization
hankL1<- function(A) {
  N<-nrow(A)+ncol(A)-1
  F<-rep(0,N)
  L.1<-min(nrow(A),ncol(A))
  K.1<-max(nrow(A),ncol(A))
  for(k in 1:N) {
    v<-A[c(row(A) + col(A) - k == 1)]
    F[k]<-median(v)
  }
  return (F)
}

#L2-hankelization
hankL2<- function(A) {
  N<-nrow(A)+ncol(A)-1
  F<-rep(0,N)
  L.1<-min(nrow(A),ncol(A))
  K.1<-max(nrow(A),ncol(A))
  for(k in 1:N) {
    v<-A[c(row(A) + col(A) - k == 1)]
    F[k]<-mean(v)
  }
  return (F)
}