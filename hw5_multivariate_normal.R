setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####4.26
#####(a)
x<-matrix(c(1:9,11,18.95,19.00,17.95,15.54,14.00,12.95,
8.94,7.49,6.00,3.99),10,2)
xbar<-matrix(colMeans(x),2,1)
s<-matrix(cov(x),2,2)
st<-solve(s)

x1<-matrix(x[1,])
x2<-matrix(x[2,])
x3<-matrix(x[3,])
x4<-matrix(x[4,])
x5<-matrix(x[5,])
x6<-matrix(x[6,])
x7<-matrix(x[7,])
x8<-matrix(x[8,])
x9<-matrix(x[9,])
x10<-matrix(x[10,])
d1<-t(x1-xbar)%*%st%*%(x1-xbar)
d2<-t(x2-xbar)%*%st%*%(x2-xbar)
d3<-t(x3-xbar)%*%st%*%(x3-xbar)
d4<-t(x4-xbar)%*%st%*%(x4-xbar)
d5<-t(x5-xbar)%*%st%*%(x5-xbar)
d6<-t(x6-xbar)%*%st%*%(x6-xbar)
d7<-t(x7-xbar)%*%st%*%(x7-xbar)
d8<-t(x8-xbar)%*%st%*%(x8-xbar)
d9<-t(x9-xbar)%*%st%*%(x9-xbar)
d10<-t(x10-xbar)%*%st%*%(x10-xbar)
d=c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)

#####(b)
length(d[d>1.39]) #得知有8個
install.packages("mixtools")
library(mixtools)
plot(x)
ellipse(mu<-colMeans(x),sigma<-cov(x),alpha=.5,col='red') #由圖形看也是

#####(c)
q=c(0.10,0.33,0.58,0.86,1.20,1.60,2.10,2.77,3.79,5.99) #自由度10卡方分布分位?
plot(q,sort(d)) #?形不是一?斜率1的直?


#####4.39
#####(a)
t406<-read.table("T4-6.DAT")
xbar1<-matrix(colMeans(t406[,1:5]),5,1)
q1<-c() #空的?列
for (i in 1:130){q1[i]<-c((i-(1/2))/130)} #造出概率水平
qnorm(q1) #在?概率水平下的分位?
par(mfcol=c(3,2))
plot(qnorm(q1),sort(t406$V1))
plot(qnorm(q1),sort(t406$V2))
plot(qnorm(q1),sort(t406$V3))
plot(qnorm(q1),sort(t406$V4))
plot(qnorm(q1),sort(t406$V5))
par(mfcol=c(1,1))

#####(b)
b<-matrix(cov(t406[,1:5]),5,5)
c<-solve(b)

t1<-array(t406$V1,c(130,1))
t2<-array(t406$V2,c(130,1))
t3<-array(t406$V3,c(130,1))
t4<-array(t406$V4,c(130,1))
t5<-array(t406$V5,c(130,1))
phy<-cbind(t1,t2,t3,t4,t5)

phy1<-array(0,c(5,130))
for (i in 1:130){phy1[,i] <- (phy[i,1:5]-xbar1)[,1]} ## difference

phy2<-array(0,c(130,1))
for (i in 1:130){phy2[i,]<-(t(phy1[,i])%*%c%*%phy1[,i])} ## distance

qchisq(q1,5) ## levels toward qi
plot(qchisq(q1,5),sort(phy2))

#####(c)
cor1<-cor(sort(phy[,1]),qnorm(q1))
cor2<-cor(sort(phy[,2]),qnorm(q1))
cor3<-cor(sort(phy[,3]),qnorm(q1))
cor4<-cor(sort(phy[,4]),qnorm(q1))
cor5<-cor(sort(phy[,5]),qnorm(q1))
cbind(cor1,cor2,cor3,cor4,cor5)

v5sq<-phy[,5]^(1/2)
cor(sort(v5sq),qnorm(q1)) ## closer to normal

#####4.40
#####(a)
t111<-read.table("T1-11.dat")
plot(t111)

(t111[,1]-array(mean(t111[,1]),c(15,1)))/var(t111[,1])
(t111[,2]-array(mean(t111[,2]),c(15,1)))/var(t111[,2]) ## find object7

#####(b)
q2<-c()
for (i in 1:15){q2[i]<-c((i-(1/2))/15)} 
v1sq<-t111[,1]^(1/2)
plot(qchisq(q2,2),sort(v1sq))

#####(c)
v2ln<-log(t111[,2])
plot(qchisq(q2,2),sort(v2ln))

#####(d)

park1<-array(v1sq,c(15,1))
park2<-array(v2ln,c(15,1))
park<-cbind(park1,park2)

xbar2<-matrix(colMeans(park[,1:2]),2,1)
par1<-array(0,c(2,15))
for (i in 1:15){par1[,i] <- (park[i,1:2]-xbar2)[,1]} ## difference

par2<-array(0,c(15,1))
for (i in 1:15){par2[i,]<-(t(par1[,i])%*%solve(cov(park))%*%par1[,i])} ## distance

plot(qchisq(q2,2),sort(par2))

###
qbar<-sum(qnorm(q1))
array(xbar1[1,],c(130,1))
###

###
r1<-(t(phy[,1]-array(xbar1[1,],c(130,1)))%*%(array(q1,c(130,1))-
array(qbar,c(130,1))))/(var(array(q1,c(130,1))-array(qbar,c(130,1)))*
var(phy[,1]-array(xbar1[1,],c(130,1))))
###

###
var(array(q1,c(130,1))-array(qbar,c(130,1))) ## var of x-xbar
var(phy[,1]-array(xbar1[1,],c(130,1))) ## var of q-qbar
###

###
a<-array(0,c(3,3))
for (i in 1:3){a[,i]<-c(2,3,1)}
###




