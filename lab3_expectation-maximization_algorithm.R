setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####exercise1
###(1)

t502 <- read.table("T5-2.dat")

install.packages("ICSNP")
library(ICSNP)

HotellingsT2(t502, mu=c(550,55,25))
## under normality assumption
## p-value is very small, so we reject

HotellingsT2(t502, mu=c(550,55,25), test="chi")
## non-normal case
## the same result

###(2)

mu <- matrix(c(550,55,25),3,1)
xmu <- colMeans(t502)
xcov <- cov(t502)

F1 <- 87*t(xmu-mu)%*%solve(xcov)%*%(xmu-mu)
1 - pf((87-3)*F1/(3*(87-1)),3,87-3) 
## p-value under normal assumption

1 - pchisq(F1,3)
## p-value under non-normal

#####exercise2
###(a)-Bon
###normal
###every variable

xbar <- matrix(c(0.766,0.508,0.438,0.161),4,1)
S <- matrix(c(0.856,0.635,0.173,0.096,0.635,0.568,
0.128,0.067,0.173,0.128,0.171,0.039,0.096,0.067,
0.039,0.043),4,4)
a <- diag(1,nrow=4)

mu1BL <- xbar[1] - qt(1-0.05/(2*4),50)*sqrt(S[1,1]/51)
mu1BU <- xbar[1] + qt(1-0.05/(2*4),50)*sqrt(S[1,1]/51)

mu2BL <- xbar[2] - qt(1-0.05/(2*4),50)*sqrt(S[2,2]/51)
mu2BU <- xbar[2] + qt(1-0.05/(2*4),50)*sqrt(S[2,2]/51)

mu3BL <- xbar[3] - qt(1-0.05/(2*4),50)*sqrt(S[3,3]/51)
mu3BU <- xbar[3] + qt(1-0.05/(2*4),50)*sqrt(S[3,3]/51)

mu4BL <- xbar[4] - qt(1-0.05/(2*4),50)*sqrt(S[4,4]/51)
mu4BU <- xbar[4] + qt(1-0.05/(2*4),50)*sqrt(S[4,4]/51)

c(mu1BL,mu1BU)
c(mu2BL,mu2BU)
c(mu3BL,mu3BU)
c(mu4BL,mu4BU)

###total

a1 <- matrix(1/4,4,1)

totalBL <- t(a1)%*%xbar - qt(1-0.05/(2*4),50)*
sqrt(t(a1)%*%S%*%a1/51)
totalBU <- t(a1)%*%xbar + qt(1-0.05/(2*4),50)*
sqrt(t(a1)%*%S%*%a1/51)
c(totalBL,totalBU)

###difference

a2 <- matrix(c(1,-1,0,0),4,1)

diffBL <- t(a2)%*%xbar - qt(1-0.05/(2*4),50)*
sqrt(t(a2)%*%S%*%a2/51)
diffBU <- t(a2)%*%xbar + qt(1-0.05/(2*4),50)*
sqrt(t(a2)%*%S%*%a2/51)
c(diffBL,diffBU)

###non-normal
###every variable

mu1NBL <- xbar[1] - qnorm(1-0.05/(2*4))*sqrt(S[1,1]/51)
mu1NBU <- xbar[1] + qnorm(1-0.05/(2*4))*sqrt(S[1,1]/51)

mu2NBL <- xbar[2] - qnorm(1-0.05/(2*4))*sqrt(S[2,2]/51)
mu2NBU <- xbar[2] + qnorm(1-0.05/(2*4))*sqrt(S[2,2]/51)

mu3NBL <- xbar[3] - qnorm(1-0.05/(2*4))*sqrt(S[3,3]/51)
mu3NBU <- xbar[3] + qnorm(1-0.05/(2*4))*sqrt(S[3,3]/51)

mu4NBL <- xbar[4] - qnorm(1-0.05/(2*4))*sqrt(S[4,4]/51)
mu4NBU <- xbar[4] + qnorm(1-0.05/(2*4))*sqrt(S[4,4]/51)

c(mu1NBL,mu1NBU)
c(mu2NBL,mu2NBU)
c(mu3NBL,mu3NBU)
c(mu4NBL,mu4NBU)

###total

totalNBL <- t(a1)%*%xbar - qnorm(1-0.05/(2*4))*
sqrt(t(a1)%*%S%*%a1/51)
totalNBU <- t(a1)%*%xbar + qnorm(1-0.05/(2*4))*
sqrt(t(a1)%*%S%*%a1/51)
c(totalNBL,totalNBU)

###difference

diffNBL <- t(a2)%*%xbar - qnorm(1-0.05/(2*4))*
sqrt(t(a2)%*%S%*%a2/51)
diffNBU <- t(a2)%*%xbar + qnorm(1-0.05/(2*4))*
sqrt(t(a2)%*%S%*%a2/51)
c(diffNBL,diffNBU)

###(b)-simul
###normal
###every variable

mu1L <- t(a[,1])%*%xbar - sqrt(S[1,1]*4*(51-1)*qf(0.95,4,47)/
(51*(51-4)))
mu1U <- t(a[,1])%*%xbar + sqrt(S[1,1]*4*(51-1)*qf(0.95,4,47)/
(51*(51-4)))

mu2L <- t(a[,2])%*%xbar - sqrt(S[2,2]*4*(51-1)*qf(0.95,4,47)/
(51*(51-4)))
mu2U <- t(a[,2])%*%xbar + sqrt(S[2,2]*4*(51-1)*qf(0.95,4,47)/
(51*(51-4)))

mu3L <- t(a[,3])%*%xbar - sqrt(S[3,3]*4*(51-1)*qf(0.95,4,47)/
(51*(51-4)))
mu3U <- t(a[,3])%*%xbar + sqrt(S[3,3]*4*(51-1)*qf(0.95,4,47)/
(51*(51-4)))

mu4L <- t(a[,4])%*%xbar - sqrt(S[4,4]*4*(51-1)*qf(0.95,4,47)/
(51*(51-4)))
mu4U <- t(a[,4])%*%xbar + sqrt(S[4,4]*4*(51-1)*qf(0.95,4,47)/
(51*(51-4)))

c(mu1L,mu1U)
c(mu2L,mu2U)
c(mu3L,mu3U)
c(mu4L,mu4U)

###total

totalL <- t(a1)%*%xbar - sqrt(t(a1)%*%S%*%a1*4*(51-1)*
qf(0.95,4,47)/(51*(51-4)))
totalU <- t(a1)%*%xbar + sqrt(t(a1)%*%S%*%a1*4*(51-1)*
qf(0.95,4,47)/(51*(51-4)))
c(totalL,totalU)

###difference

diffL <- t(a2)%*%xbar - sqrt(t(a2)%*%S%*%a2*4*(51-1)*
qf(0.95,4,47)/(51*(51-4)))
diffU <- t(a2)%*%xbar + sqrt(t(a2)%*%S%*%a2*4*(51-1)*
qf(0.95,4,47)/(51*(51-4)))
c(diffL,diffU)

###non-normal
###every variable

mu1NL <- xbar[1] - sqrt(qchisq(0.95,4)*S[1,1]/51)
mu1NU <- xbar[1] + sqrt(qchisq(0.95,4)*S[1,1]/51)

mu2NL <- xbar[2] - sqrt(qchisq(0.95,4)*S[2,2]/51)
mu2NU <- xbar[2] + sqrt(qchisq(0.95,4)*S[2,2]/51)

mu3NL <- xbar[3] - sqrt(qchisq(0.95,4)*S[3,3]/51)
mu3NU <- xbar[3] + sqrt(qchisq(0.95,4)*S[3,3]/51)

mu4NL <- xbar[4] - sqrt(qchisq(0.95,4)*S[4,4]/51)
mu4NU <- xbar[4] + sqrt(qchisq(0.95,4)*S[4,4]/51)

c(mu1NL,mu1NU)
c(mu2NL,mu2NU)
c(mu3NL,mu3NU)
c(mu4NL,mu4NU)

###total

totalNL <- t(a1)%*%xbar - sqrt(t(a1)%*%S%*%a1*qchisq(0.95,4)/51)
totalNU <- t(a1)%*%xbar + sqrt(t(a1)%*%S%*%a1*qchisq(0.95,4)/51)
c(totalNL,totalNU)

###difference

diffNL <- t(a2)%*%xbar - sqrt(t(a2)%*%S%*%a2*qchisq(0.95,4)/51)
diffNU <- t(a2)%*%xbar + sqrt(t(a2)%*%S%*%a2*qchisq(0.95,4)/51)
c(diffNL,diffNU)

#####exercise3

em<-function(xdata,mu0,Sigma0){
  n<-nrow(xdata)
  p<-ncol(xdata)

   err<-function(mu0,Sigma0,mu1,Sigma1){
    th0<-c(mu0,as.vector(Sigma0))
    th1<-c(mu1,as.vector(Sigma1))
    sqrt(sum((th0-th1)*(th0-th1)))
   }

  mu1<-mu0+1
  Sigma1<-Sigma0+1
  while(err(mu0,Sigma0,mu1,Sigma1)>1e-6){
     mu1<-mu0
     Sigma1<-Sigma0

     zdata<-xdata
     Ai<-matrix(0,p,p)
    for(i in 1:n){
     if(any(is.na(xdata[i,]))){
      zi<-xdata[i,]
      na.idx<-(1:p)[is.na(zi)]
      cs.idx<-(1:p)[-na.idx]
      Sigma012<-Sigma0[na.idx,cs.idx,drop=FALSE]
      Sigma022.iv<-solve(Sigma0[cs.idx,cs.idx])
      zdata[i,na.idx]<-mu0[na.idx]+(Sigma012%*%Sigma022.iv)%*%(zi[cs.idx]-mu0[cs.idx])
      Ai[na.idx,na.idx]<-Ai[na.idx,na.idx]+Sigma0[na.idx,na.idx]-Sigma012%*%Sigma022.iv%*%t(Sigma012)
       }
     }
    mu0<-colMeans(zdata)
    Sigma0<-(n-1)*cov(zdata)/n+Ai/n
   }
 return(list(mu=mu0,Sigma=Sigma0))
}

csvtao <- read.csv("tao.csv", row.names=1)

library(MASS)

mtao <- as.matrix(csvtao)

mu0 <- apply(mtao[,4:8],2,mean,na.rm=TRUE)
nas <- is.na(mtao[,4:8])
na.num <- apply(nas,2,sum)
zdata <- mtao[,4:8]
zdata[nas] <- rep(mu0,na.num)
Sigma0 <- cov(zdata)

system.time(rlt <- em(mtao[,4:8],mu0,Sigma0))
rlt

library(xtable)

rlt.mu <- as.matrix(rlt$mu)
rlt.Sigma <- as.matrix(rlt$Sigma)

xtable(rlt.mu)
xtable(rlt.Sigma)

