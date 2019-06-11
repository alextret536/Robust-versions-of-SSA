library(Rssa)
library(pcaL1)
library(pcaMethods)
library(matrixcalc)
#N <- 500
N<-240
Per<-120
L<-120
sig <-(1:N)*exp(4*(1:N)/(N))*sin(2*pi*(1:N)/30)

plot(sig,type='l')
X<-hankel(sig,L=120)

M <- 10

MSE.SSA <- rep(0,M)
MSE.L1svd <- rep(0,M)
MSE.IRLS <- rep(0,M)
MSE.IRLSmodif <- rep(0,M)
#MSE.series <- rep(0,M)

MAD.SSA <- rep(0,M)
MAD.L1svd <- rep(0,M)
MAD.IRLS <- rep(0,M)
MAD.IRLSmodif <- rep(0,M)
#MAD.series <- rep(0,M)

timeSSA<-0
timePCAL1<-0
timeRodrigues<-0
timeIRLS<-0

for(k in 1:M)
{
  
  sig <-(1:N)*exp(4*(1:N)/(N))*sin(2*pi*(1:N)/30)
  sig.outl<-sig
  outlier.seq<-sample(1:N,N*0.01)
  sig.outl[outlier.seq]<-sig.outl[outlier.seq]+1.5*sig.outl[outlier.seq]
  #sig.outl[outlier.seq]<-sig.outl[outlier.seq]+1.5*sig.outl[outlier.seq]
  ser <- sig.outl + rnorm(N)
  plot(ser,type='l')
  title(main="Series plot")
  
  
  ############################### #SSA120.3 ##############################
  start.time<-Sys.time()
  s <- ssa(ser, L = 120)
  rec <- reconstruct(s, groups = list(c(1:4))) 
  end.time<-Sys.time()
  timeSSA<-timeSSA+(end.time-start.time)
  # plot(s, type = "vectors",idx = 1:10)
  # plot(rec, add.residuals = TRUE, add.original = TRUE,
  #  plot.method = "xyplot", superpose = TRUE, auto.key = list(columns = 2))
  trend.season <- rec$F1
  
  ############################# #pcaL1(120.3)###############################
  X<-hankel(ser,L=120)
  start.time<-Sys.time()
  s120.L1svd<-l1pca(X,center=FALSE,projections="l1",projDim=4)
  end.time<-Sys.time()
  Pr<-s120.L1svd$projPoints
  Pr120.L1svd<-hankL1(Pr)
  timePCAL1<-timePCAL1+(end.time-start.time)
  
  ###################### WLS ########################
  start.time<-Sys.time()
  Pr<-QR_WLS(X,4)
  end.time<-Sys.time()
  Pr120.IRLS<-hankL1(Pr)
  timeIRLS<-timeIRLS+(end.time-start.time)
  
  ###################### WLS (modif.) ########################
  Pr<-QR_WLS1(X,4)
  Pr120.IRLSmodif<-hankL1(Pr)
  
  
  # ############################ #rodr ######################################
  # start.time<-Sys.time()
  # s<-robustSvd(X)
  # d.decr<-s$d[order(s$d,decreasing = TRUE)]
  # u.decr<-s$u[,order(s$d,decreasing=TRUE)]
  # v.decr<-s$v[,order(s$d,decreasing=TRUE)]
  # Pr<-d.decr[1]*u.decr[,1]%*%t(v.decr[,1])+d.decr[2]*u.decr[,2]%*%t(v.decr[,2])+ d.decr[3]*u.decr[,3]%*%t(v.decr[,3])+d.decr[4]*u.decr[,4]%*%t(v.decr[,4])
  # end.time<-Sys.time()
  # Pr120.L1svd.rodrigues<-hankL1(Pr)
  # timeRodrigues<-timeRodrigues+(end.time-start.time)
  
  #MSE 
  MSE.SSA[k] <- mean((sig - trend.season)[1:N]^2)
  MSE.L1svd[k] <- mean((sig - Pr120.L1svd)[1:N]^2)
  MSE.IRLS[k] <- mean((sig - Pr120.IRLS)[1:N]^2)
  MSE.IRLSmodif[k] <- mean((sig - Pr120.IRLSmodif)[1:N]^2)
  #MSE.series[k] <- mean((sig - ser)[1:N]^2)
  
  #MAD
  MAD.SSA[k] <- mean(abs((sig - trend.season)[1:N]))
  MAD.L1svd[k] <- mean(abs((sig - Pr120.L1svd)[1:N]))
  MAD.IRLS[k] <- mean(abs((sig - Pr120.IRLS)[1:N]))
  MAD.IRLSmodif[k] <- mean(abs((sig - Pr120.IRLSmodif)[1:N]))
  #MAD.series[k] <- mean(abs((sig - ser)[1:N]))
  
}


RMSE.SSA<-sqrt(mean(MSE.SSA))
RMSE.L1svd<-sqrt(mean(MSE.L1svd)) 
RMSE.IRLS<-sqrt(mean(MSE.IRLS))
RMSE.IRLSmodif<-sqrt(mean(MSE.IRLSmodif))
#RMSESSA120.series<-sqrt(mean(MSE.series))

MAD.SSA<-mean(MAD.SSA)
MAD.L1svd<-mean(MAD.L1svd)
MAD.IRLS<-mean(MAD.IRLS)
MAD.IRLSmodif<-mean(MAD.IRLSmodif)
#MADSSA120.series<-mean(MAD.series)

