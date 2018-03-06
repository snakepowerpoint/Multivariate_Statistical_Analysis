setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####5.4
###(a)

t504<-read.table("T5-1.dat")
sigma<-cov(t504)
en=eigen(sigma)
en$vectors ## axes of ellipsoid

F504<-qf(p=0.9, df1=3, df2=17)
F504a<-F504*(3*(20-1)/(20-3))

sqrt(en$values[1])*sqrt(F504a/20) ## lengths of axis1
sqrt(en$values[2])*sqrt(F504a/20) ## lengths of axis2
sqrt(en$values[3])*sqrt(F504a/20) ## lengths of axis3

###(b)

qqnorm(t504$V1) ## check normal assumption
qqnorm(t504$V2)
qqnorm(t504$V3)

plot(t504) ## plot pairs of observations

#####5.7
###simultaneous T^2 intervals

a1=matrix(c(1,0,0),3,1)
a2=matrix(c(0,1,0),3,1)
a3=matrix(c(0,0,1),3,1)

xbar<- apply(t504,2,mean) ## or xbar<- colMeans(t504)
F507<-qf(p=0.95, df1=3, df2=17)
F507a<-F507*(3*(20-1)/(20-3))

low1<- t(a1)%*%xbar - sqrt((t(a1)%*%sigma%*%a1)*F507a/20)
upp1<- t(a1)%*%xbar + sqrt((t(a1)%*%sigma%*%a1)*F507a/20)
range1<- c(low1,upp1)

low2<- t(a2)%*%xbar - sqrt((t(a2)%*%sigma%*%a2)*F507a/20)
upp2<- t(a2)%*%xbar + sqrt((t(a2)%*%sigma%*%a2)*F507a/20)
range2<- c(low2,upp2)

low3<- t(a3)%*%xbar - sqrt((t(a3)%*%sigma%*%a3)*F507a/20)
upp3<- t(a3)%*%xbar + sqrt((t(a3)%*%sigma%*%a3)*F507a/20)
range3<- c(low3,upp3)

###Bonferroni intervals

T507<- qt(p=1-0.05/(2*3), df=19)
var<- array(diag(sigma))

low11<- xbar[1] - sqrt(var[1]/20)*T507
upp11<- xbar[1] + sqrt(var[1]/20)*T507
range11<- c(low11,upp11)

low22<- xbar[2] - sqrt(var[2]/20)*T507
upp22<- xbar[2] + sqrt(var[2]/20)*T507
range22<- c(low22,upp22)

low33<- xbar[3] - sqrt(var[3]/20)*T507
upp33<- xbar[3] + sqrt(var[3]/20)*T507
range33<- c(low33,upp33)

## compare simultaneous T^2 intervals with Bonferronis',
## we find the Bonferroni intervals are more narrower

#####5.12
## the function follower is used to run EM algorithm

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

em512<-matrix(c(3,4,NA,5,6,4,8,NA,0,3,3,NA),4,3) 
mu0<-apply(em512,2,mean,na.rm=TRUE) ## the initial guess of mean
nas<-is.na(em512)
na.num<-apply(nas,2,sum)
zdata<-em512
zdata[nas]<-rep(mu0,na.num) ## replace NA with the initial mean
Sigma0<-(4-1)*cov(zdata)/4 ## the initial guess of variance(n)
rlt<-em(em512,mu0,Sigma0) ## the result of EM algorithm

system.time(rlt<-em(em512,mu0,Sigma0)) ## how much time have been used

#####5.20
###(a)

install.packages("mixtools")
library(mixtools)

t520<- read.table("T5-12.dat")

conf.reg<-function(xdata,alpha){
  if(ncol(xdata)!=2) stop("Only for bivariate normal")
  n<-nrow(xdata)
  xbar<-colMeans(xdata)
  S<-cov(xdata)
  es<-eigen(S)
  e1<-es$vec %*% diag(sqrt(es$val))
  r1<-sqrt(qf(alpha,2,n-2))*sqrt(2*(n-1)/(n*(n-2)))
  theta<-seq(0,2*pi,len=250)
  v1<-cbind(r1*cos(theta), r1*sin(theta))
  pts<-t(xbar-(e1%*%t(v1)))
  plot(pts,type="l",main="Confidence Region for Bivariate Normal",xlab=colnames(xdata)[1],ylab=colnames(xdata)[2],asp=1)
  segments(0,xbar[2],xbar[1],xbar[2],lty=2) # highlight the center
  segments(xbar[1],0,xbar[1],xbar[2],lty=2)
  
  th2<-c(0,pi/2,pi,3*pi/2,2*pi)   #adding the axis
  v2<-cbind(r1*cos(th2), r1*sin(th2))
  pts2<-t(xbar-(e1%*%t(v2)))
  segments(pts2[3,1],pts2[3,2],pts2[1,1],pts2[1,2],lty=3)  
  segments(pts2[2,1],pts2[2,2],pts2[4,1],pts2[4,2],lty=3)
}  ## ellipsoid function

conf.reg(t520, .95)
points(t(apply(t520,2,mean)),col='red',pch=19) ## the center of ellipse
points(190,275,pch=20,col='red') ## the point is in the ellipse

## (190,275) is in the ellipse

###(b)
###simultaneous T^2 intervals

c1=matrix(c(1,0),2,1)
c2=matrix(c(0,1),2,1)
xbar1<- colMeans(t520)
F520<- qf(p=0.95, df1=2, df2=43)
F520a<- F520*(3*(45-1)/(45-3))
sigma1<- cov(t520) 

lowa<- t(c1)%*%xbar1 - sqrt((t(c1)%*%sigma1%*%c1)*F520a/45)
uppa<- t(c1)%*%xbar1 + sqrt((t(c1)%*%sigma1%*%c1)*F520a/45)
rangea<- c(lowa,uppa)

lowb<- t(c2)%*%xbar1 - sqrt((t(c2)%*%sigma1%*%c2)*F520a/45)
uppb<- t(c2)%*%xbar1 + sqrt((t(c2)%*%sigma1%*%c2)*F520a/45)
rangeb<- c(lowb,uppb)

###Bonferroni intervals

T520<- qt(p=1-0.05/(2*2), df=44)
var1<- array(diag(sigma1))

lowaa<- xbar1[1] - sqrt(var1[1]/45)*T520
uppaa<- xbar1[1] + sqrt(var1[1]/45)*T520
rangeaa<- c(lowaa,uppaa)

lowbb<- xbar1[2] - sqrt(var1[2]/45)*T520
uppbb<- xbar1[2] + sqrt(var1[2]/45)*T520
rangebb<- c(lowbb,uppbb)

## compare simultaneous T^2 intervals with Bonferronis',
## we find the Bonferroni intervals are more narrower

###(c)

qqnorm(t520$V1) ## check normal assumption
qqnorm(t520$V2) 
plot(t520$V1,t520$V2)
plot(t520$V1,t520$V2,xlim=c(160,230),ylim=c(220,320))

D1<-mahalanobis(t520,colMeans(t520), cov(t520)) ## square of distance and sort

qqplot(qchisq(ppoints(45), df = 2), D1,
main = expression("Q-Q plot for" ~~ {chi^2}[nu == 2]))
abline(c(0,1), col='red')

## we can claim that the normal assumption is reasonable



