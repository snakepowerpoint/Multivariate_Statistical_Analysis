setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####1.6
###(a)

t105 <- read.table("T1-5.dat")
names(t105)
names(t105) <- c("wind","solar","CO","NO","NO2","O3","HC")

par(mfrow=c(4,2))
stripchart(t105[,1], "stack", xlab="wind")
stripchart(t105[,2], "stack", xlab="solar")
stripchart(t105[,3], "stack", xlab="CO")
stripchart(t105[,4], "stack", xlab="NO")
stripchart(t105[,5], "stack", xlab="NO2")
stripchart(t105[,6], "stack", xlab="O3")
stripchart(t105[,7], "stack", xlab="HC")

###(b)

library(xtable)

colMeans(t105)
cov(t105)
cor(t105)
## CO,NO,NO2,O3 have large correlation

xtable(as.matrix(colMeans(t105)))
xtable(cov(t105))
xtable(cor(t105))


#####1.9
###(a)

x1 <- c(-6,-3,-2,1,2,5,6,8)
x2 <- c(-2,-3,1,-1,2,1,5,3)
plot(x1,x2)

var(x1)
var(x2)
cov(x1,x2)

###(b)

r <- matrix(c(0.899,-0.438,0.438,0.899),2,2)
x <- rbind(x1,x2)	
(r%*%x)[1,]
(r%*%x)[2,]

xtable(t(as.matrix((r%*%x)[1,])))
xtable(t(as.matrix((r%*%x)[2,])))

###(c)

x11 <- (r%*%x)[1,]
x22 <- (r%*%x)[2,]
var(x11)
var(x22)

###(d)

p <- r%*%matrix(c(4,-2),2,1)
sqrt(p[1,]^2/var(x11)+p[2,]^2/var(x22))

###(e)

var(x1)*0.899^2 + 2*0.899*0.438*cov(x1,x2) + var(x2)*0.438^2
var(x1)*0.438^2 - 2*0.899*0.438*cov(x1,x2) + var(x2)*0.899^2

a11 <- 0.899^2/(var(x1)*0.899^2 + 2*0.899*0.438*cov(x1,x2) + var(x2)*0.438^2)+
0.438^2/(var(x1)*0.438^2 - 2*0.899*0.438*cov(x1,x2) + var(x2)*0.899^2)
a22 <- 0.438^2/(var(x1)*0.899^2 + 2*0.899*0.438*cov(x1,x2) + var(x2)*0.438^2)+
0.899^2/(var(x1)*0.438^2 - 2*0.899*0.438*cov(x1,x2) + var(x2)*0.899^2)
a12 <- 0.438*0.899/(var(x1)*0.899^2 + 2*0.899*0.438*cov(x1,x2) + var(x2)*0.438^2)-
0.438*0.899/(var(x1)*0.438^2 - 2*0.899*0.438*cov(x1,x2) + var(x2)*0.899^2)
sqrt(a11*4^2+2*a12*4*(-2)+a22*(-2)^2) 
## the same


#####1.13

x<-c("E","D","C","B","A")
y<-c(1:5)
n<-1000
p1 <- 0;
p2 <- 0;
for(i in 1:length(x)){
for(j in 1:length(y)){
m<-abs(i-5)+abs(i-3)+abs(i-1)+abs(j-2)+abs(j-5)+abs(j-3)
if(m<n){p1<-i;p2<-j;n<-m}
}}
x[p1];y[p2]


#####3.9
###(a)

X <- matrix(c(12,18,14,20,16,17,20,16,18,19,29,38,30,38,35),5,3)
Xbar <- colMeans(X)
X.1 <- X - t(matrix(rep(Xbar,5),3,5))
## a'=(1,1,-1)

###(b)

xtable(cov(X))
solve(cov(X)) ## 0
cov(X)%*%matrix(c(1,1,-1),3,1) ## 0

###(c)
## as above


#####3.14
###(a)

X <- matrix(c(9,5,1,1,3,2),3,2)
c <- matrix(c(-1,2),2,1)
b <- matrix(c(2,3),2,1)

cX <- t(c)%*%t(X)
bX <- t(b)%*%t(X)

mean(cX)
mean(bX)
var(as.vector(cX))
var(as.vector(bX))
cov(as.vector(cX), as.vector(bX))

###(b)

colMeans(X)%*%c
colMeans(X)%*%b

t(c)%*%cov(X)%*%c
t(b)%*%cov(X)%*%b
t(c)%*%cov(X)%*%b
## the same


#####3.20
###(a)

t302 <- read.table("T3-2.dat")
names(t302) <- c("SnowTime", "LaborTime")

mean(t302[,1]-t302[,2])
var(t302[,1]-t302[,2])

###(b)

mean(t302[,2]-t302[,1]) ## 
var(t302[,2]-t302[,1]) ## the same







