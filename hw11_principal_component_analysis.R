setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####8.6
###(a)

S.86 <- matrix(c(7476.45,303.62,303.62,26.19),2,2)
eigen(S.86)

###(b)

eigen(S.86)$values[1]/(eigen(S.86)$values[1]+eigen(S.86)$values[2])

###(c)

library(mixtools)

x1 <- c(108.28,152.36,95.04,65.45,62.97,263.99,265.19,285.06,92.01,165.68)
x2 <- c(17.05,16.59,10.91,14.14,9.52,25.33,18.54,15.73,8.10,11.13)
x3 <- c(1484.10,750.33,766.42,1110.46,1031.29,195.26,193.83,191.11,1175.16,211.15)
t104 <- cbind(x1,x2,x3)
t104 <- data.frame(t104)
names(t104)=c("Sales","Profits","Assets")

mu.86 <- c(155.60,14.70)
pchisq(1.4,2)

plot(t104[,1:2], type="n", xlim=c(0,300), main="50% contour line")
ellipse(mu=mu.86, sigma=S.86, alpha=0.5)
points(t(mu.86), col='red', pch=19)
points(t104[,1:2])
sqrt(eigen(S.86)$values[1]*1.4) ## half length of axis1
sqrt(eigen(S.86)$values[2]*1.4) ## half length of axis2

###(d)

eigen(S.86)$vec[1,1]*sqrt(eigen(S.86)$val[1])/sqrt(S.86[1,1])
eigen(S.86)$vec[1,2]*sqrt(eigen(S.86)$val[1])/sqrt(S.86[2,2])
## the first component is almost determined by x1
## we can get the same conclusion by correlation between x1 and y1

#--------------------------------test------------------------------

pc.104 <- princomp(t104)
summary(pc.104, loadings=T)
eigen(cov(t104))
attributes(pc.104)

test1 <- t(matrix(pc.104$loadings[,1]))
test2 <- matrix(as.matrix(t104[1,]))
test1
test2
test2 <- test2 - colMeans(t104)
test1%*%test2
pc.104$scores ## exactly the same
## so we know that "scores" in R is (ei)'%*%(Xj - Xbar)

xbar <- colMeans(t104)
xbar <- matrix(as.matrix(xbar),1,3)
t104.bar <- t104 - matrix(rbind(xbar),nrow=10,ncol=3,byrow=T)
t104.bar
cov(t104)
cov(t104.bar) ## the variance of (x-xbar) is same as x

#--------------------------------test------------------------------

#####8.10
###(a)

t804 <- read.table("T8-4.dat", col.names=c("JPMorgan","CiTibank",
"WellsFargo","RoyalDutchShell","ExxonMobil"))

S.804 <- cov(t804)
m.804 <- colMeans(t804)
pc.804 <- princomp(t804)
summary(pc.804, loadings=T)

library(xtable)
xtable(S.804) ## too small
loading.804 <- matrix(loadings(pc.804),5,5)
xtable(loading.804)

###(b)

summary(pc.804, loadings=T)
## component1 represents market component
## component2 represents industry component
## component3 is difficult to explain

###(c)

evalue.CI <- function(datamatrix, labelvec, m=length(labelvec), conf.level=0.95){
alpha<- 1-conf.level
n<-nrow(datamatrix)
z<-qnorm(alpha/(2*m),lower=F)
lambdas<- princomp(datamatrix)$sdev^2
print(lambdas)
LCL<-lambdas[labelvec]/(1+z*sqrt(2/n))
UCL<-lambdas[labelvec]/(1-z*sqrt(2/n))
CIs<-cbind(LCL,UCL)
return(CIs)
}

evalue.CI(t804,c(1,2,3),conf.level=0.9)

###(d)

S.804
eigen(S.804)
## eigen values are not all the same, and non-diagonal elements
## are not all zero. So we can express this data under
## smaller dimension

#####8.18
###(a)

t109 <- read.table("T1-9.dat")
names(t109)
names(t109) <- c("Country","100m","200m","400m","800m.min","1500m.min",
"3000m.min","marathon.min")

R.109 <- cor(t109[,2:8])
eigen(R.109)

###(b)

t109.s <- scale(t109[,2:8])
pc.109.s <- princomp(t109.s)
summary(pc.109.s, loadings=T) ## cumulative proportion 91.95%
cor(t109.s)
R.109 ## the same

###(c)

summary(pc.109.s, loadings=T)

###(d)

score.109 <- data.frame(t109$Country,pc.109.s$scores[,1])
names(score.109) <- c("Country","Scores")
score.109[order(score.109[,2],decreasing=T),]

#####8.19


#####8.20
###(a)

t806 <- read.table("T8-6.dat")
names(t806)
names(t806) <- c("Country","100m","200m","400m","800m.min","1500m.min",
"5000m.min","10000m.min","marathon.min")

R.806 <- cor(t806[,2:9])
eigen(R.806)

xtable(t(eigen(R.806)$vectors))
value.806 <- matrix(eigen(R.806)$values, 1,8)
xtable(value.806)

###(b)

t806.s <- scale(t806[,2:9])
pc.806.s <- princomp(t806.s)
summary(pc.806.s, loadings=T) ## cumulative proportion 91.77%
cor(t806.s)
R.806 ## the same

xtable(matrix(pc.806.s$loadings,8,8))
xtable(cor(t806.s))

###(c)

summary(pc.806.s, loadings=T)
## component1 represents every type of contest
## component2 represents difference between short and long distance

###(d)

score.806 <- data.frame(t806$Country,pc.806.s$scores[,1])
names(score.806) <- c("Country","Scores")
contrast <- data.frame(score.806[order(score.806[,2],decreasing=T),],
score.109[order(score.109[,2],decreasing=T),])
contrast
## we find that there is some difference between men and women

xtable(contrast[1:10,])

#####8.22

###using S
###(a)

t110 <- read.table("T1-10.dat")
names(t110)
names(t110) <- c("Breed","SalePr","YrHgt","FtFrBody","PrctFFB",
"Frame","BkFat","SaleHt","SaleWt")

head(t110)
pc.110 <- princomp(t110[,3:9])
summary(pc.110, loadings=T)
cov(t110)[,3:9]

plot(1:(length(pc.110$sdev)), (pc.110$sdev)^2, type='b', 
main="Scree Plot", xlab="Number of Components", ylab="Eigenvalue Size")
## choose the first two components
screeplot(pc.110) ## the similar method

###(b)

summary(pc.110, loadings=T)
## component1 represents FtFrBody & SaleWt
## component2 represents difference between FtFrBody and SaleWt

loading.110 <- matrix(pc.110$loadings,7,7)
xtable(loading.110[,1:2])

###(c)

## we can get such index because the body size and configuration of ox
## are almost equivalent of FtFrBody and SaleWt

###(d)

plot(pc.110$scores[,1:2], type="n")
text(pc.110$scores[,1], pc.110$scores[,2], labels=c(1:76), cex=0.7, lwd=2,
col=c(rep("red", times=32), rep("blue", times=17), rep("green", times=27)))
## outliers 50,51

###(e)

attributes(pc.110)
qqnorm(pc.110$scores[,1]) ## nearly normal distribution
qqline(pc.110$scores[,1])

###using R
###(a)

pc.110.r <- princomp(t110[,3:9], cor=T)
summary(pc.110.r, loadings=T)

plot(1:(length(pc.110.r$sdev)), (pc.110.r$sdev)^2, type='b', 
main="Scree Plot", xlab="Number of Components", ylab="Eigenvalue Size")
## choose the first two components
screeplot(pc.110.r) ## the similar method

###(b)

summary(pc.110.r, loadings=T)
## component1 represents average except for BkFat
## component2 represents difference between PrctFFB and BkFat & SaleWt
## component3 represents difference between YrHgt & Frame and
## FtFrBody & PrctFFB

loading.110.r <- matrix(pc.110.r$loadings,7,7)
xtable(loading.110.r[,1:2])

###(c)

## we can get such index because the body size and configuration of ox
## can be defined by principal component

###(d)

plot(pc.110.r$scores[,1:2], type="n")
text(pc.110.r$scores[,1], pc.110.r$scores[,2], labels=c(1:76), cex=0.7,
lwd=2, col=c(rep("red", times=32), rep("blue", times=17),
rep("green", times=27)))
## outliers 16,51

###(e)

qqnorm(pc.110.r$scores[,1]) ## nearly normal distribution
qqline(pc.110.r$scores[,1])

#------------------------------------------------------------
john <- cov(t110[,3:9])
eigen(john)
#------------------------------------------------------------
## scores VS y ??
#------------------------------------------------------------

summary(pc.110, loadings=T)
john1 <- as.matrix(pc.110$loadings[,1])
john2 <- as.matrix(t110[1,3:9])
a <- matrix(john1)
b <- matrix(john2)

t(a)%*%b
pc.110$scores[1,1] ## why are they different?

#-----------------------test---------------------------------
pc.86 <- princomp(t104[,1:2])
summary(pc.86, loadings=T)

attributes(pc.86)
pc.86$scores

en1 <- eigen(S.86)$vectors[1,]
en2 <- eigen(S.86)$vectors[2,]
en1 <- as.matrix(en1)
en2 <- as.matrix(en2)
m104 <- as.matrix(t104)

t(en1)%*%m104[1,1:2]

pc.86

e1 <- eigen(S.86)$vec %*% diag(sqrt(eigen(S.86)$val))
r1 <- sqrt(qchisq(0.5,2))
th2 <- c(0,pi/2,pi,3*pi/2,2*pi)   #adding the axis
v2 <- cbind(r1*cos(th2), r1*sin(th2))
pts2 <- t(mu.86-(e1%*%t(v2)))
segments(pts2[3,1],pts2[3,2],pts2[1,1],pts2[1,2],lty=3)  
segments(pts2[2,1],pts2[2,2],pts2[4,1],pts2[4,2],lty=3)


