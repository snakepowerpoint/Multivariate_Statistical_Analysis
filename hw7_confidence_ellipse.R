setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####6.11
## do with pen and paper

#####6.12

n1.612<- 30
n2.612<- 30
xbar1.612<- matrix(c(6.4,6.8,7.3,7.0),4,1)
xbar2.612<- matrix(c(4.3,4.9,5.3,5.1),4,1)
sp.612<- matrix(c(0.61,0.26,0.07,0.16,0.26,0.64,0.17,0.14,0.07,
0.17,0.81,0.03,0.16,0.14,0.03,0.31),4,4)
c.612<- matrix(c(1,0,-2,1,1,-2,0,1),2,4)

t((xbar1.612 + xbar2.612))%*%t(c.612)%*%solve(c.612%*%sp.612
%*%t(c.612))%*%c.612%*%(xbar1.612 + xbar2.612)*(1/n1.612 + 1/n2.612)

qf(0.95,2,57) ## can not reject

#####6.19
###(a)

t610<- read.table("T6-10.dat")
gas<- rbind(t610[1:36,1:3])
die<- rbind(t610[37:59,1:3])
cov(gas)
cov(die) ## var-cov matrix is different

mugas<- colMeans(gas)
mudie<- colMeans(die)

t(mugas-mudie)%*%solve(1/36*cov(gas)+1/23*cov(die))%*%(mugas-mudie)
qchisq(0.99,3) ## reject hypothesis

###(b)

solve(1/36*cov(gas)+1/23*cov(die))%*%(mugas-mudie) ## coefficient of a

###(c)

c.619<- matrix(c(1,0,0,0,1,0,0,0,1),3,3)
c.619[1,]%*%(mugas-mudie)-sqrt(qchisq(0.99,3)*
c.619[1,]%*%(1/36*cov(gas)+1/23*cov(die))%*%c.619[1,]) 
## mu1 lower
c.619[1,]%*%(mugas-mudie)+sqrt(qchisq(0.99,3)*
c.619[1,]%*%(1/36*cov(gas)+1/23*cov(die))%*%c.619[1,]) 
## mu1 upper

c.619[2,]%*%(mugas-mudie)-sqrt(qchisq(0.99,3)*
c.619[2,]%*%(1/36*cov(gas)+1/23*cov(die))%*%c.619[2,]) 
## mu2 lower
c.619[2,]%*%(mugas-mudie)+sqrt(qchisq(0.99,3)*
c.619[2,]%*%(1/36*cov(gas)+1/23*cov(die))%*%c.619[2,]) 
## mu2 upper

c.619[3,]%*%(mugas-mudie)-sqrt(qchisq(0.99,3)*
c.619[3,]%*%(1/36*cov(gas)+1/23*cov(die))%*%c.619[3,]) 
## mu3 lower
c.619[3,]%*%(mugas-mudie)+sqrt(qchisq(0.99,3)*
c.619[3,]%*%(1/36*cov(gas)+1/23*cov(die))%*%c.619[3,]) 
## mu3 upper
## we find variable3 is different

###(d)
install.packages("mvoutlier")
library(mvoutlier)

aq.plot(gas) ## should remove sample 9 and 21

rgas=gas[-c(9,21),]
murgas<- colMeans(rgas)

t(murgas-mudie)%*%solve(1/34*cov(rgas)+1/23*cov(die))%*%(murgas-mudie)
qchisq(0.99,3) ## reject hypothesis too

#####6.20
###(a)

t611<- read.table("T6-11.dat")
t512<- read.table("T5-12.dat")

plot(t611)
aq.plot(t611) ## remove sample 31

###(b)
## remove it 

rt611=t611[-31,]
cov(rt611)
cov(t512) ## var-cov matrix is different

murmale<- colMeans(rt611)
mufe<- colMeans(t512)

t(murmale-mufe)%*%solve(1/44*cov(rt611)+1/45*cov(t512))%*%(murmale-mufe)
qchisq(0.95,2) ## reject hypothesis

solve(1/44*cov(rt611)+1/45*cov(t512))%*%(murmale-mufe) ## coefficient of a

## modify it

mt611<- t611
fix(mt611)
cov(mt611)
cov(t512) ## var-cov matrix is different

mummale<- colMeans(mt611)
t(mummale-mufe)%*%solve(1/45*cov(mt611)+1/45*cov(t512))%*%(mummale-mufe)
qchisq(0.95,2) ## reject hypothesis

solve(1/44*cov(mt611)+1/45*cov(t512))%*%(mummale-mufe) ## similar to above

###(c)
## use modified sample

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
}

diff<- mt611-t512
colnames(diff)=c("tail","wing")
conf.reg(diff,0.95)

c.620<- matrix(c(1,0,0,1),2,2)
c.620[1,]%*%(mummale-mufe)-sqrt(qchisq(0.95,2)*
c.620[1,]%*%(1/45*cov(mt611)+1/45*cov(t512))%*%c.620[1,]) 
## mu1 lower
c.620[1,]%*%(mummale-mufe)+sqrt(qchisq(0.95,2)*
c.620[1,]%*%(1/45*cov(mt611)+1/45*cov(t512))%*%c.620[1,]) 
## mu1 upper

c.620[2,]%*%(mummale-mufe)-sqrt(qchisq(0.95,2)*
c.620[2,]%*%(1/45*cov(mt611)+1/45*cov(t512))%*%c.620[2,]) 
## mu2 lower
c.620[2,]%*%(mummale-mufe)+sqrt(qchisq(0.95,2)*
c.620[2,]%*%(1/45*cov(mt611)+1/45*cov(t512))%*%c.620[2,]) 
## mu2 upper

###(d)
## tails of female are larger than male


##############
library(mixtools)

plot(diff,type="n",xlim=c(-15,0),ylim=c(-12,12))
ellipse(colMeans(diff),1/45*cov(mt611)+1/45*cov(t512),
alpha=0.05,col='red')
## another method to plot confidence ellipse


