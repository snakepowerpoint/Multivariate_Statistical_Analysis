setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####exercise1
install.packages("mixtools")
library(mixtools)

mu<- c(1,0,0,1,1)
rho<- toeplitz(c(1,0.5,(0.5)^2,(0.5)^3,(0.5)^4))
n<-1000

library(MASS)

biv<- mvrnorm(n,mu,rho)

install.packages("car")
library(car)

pairs(biv) ## plot all variables

scatterplotMatrix(biv, ellipse=TRUE, levels=c(0, .95)) ## 95% contour
spm(biv) ## 95% contour

######exercise2

install.packages("energy")
library(energy)
install.packages("mvoutlier")
library(mvoutlier)
install.packages("rgl")
library(rgl)

t110<- read.table("T1-10.dat")
colnames(t110)<- c("breed","salepr","yrhgt","ftfrbody","prctffb",
"frame","bkfat","saleht","salewt")

t110a<- cbind(t110[,3],t110[,4],t110[,5],t110[,7],t110[,8],t110[,9])

pairs(t110a) ## check normal assumption
D0<-mahalanobis(t110a,colMeans(t110a),cov(t110a))
qqplot(qchisq(ppoints(76), df = 5), D0,
main = expression("Q-Q plot for" ~~ {chi^2}[nu == 2]))
abline(c(0,1)) ## chi-squre qq plot

aq.plot(t110a) ## find many samples....
chisq.plot(t110a) ## fine sample 16 and 51

######exercise3
###(1)

t302<- read.table("T3-2.dat")
plot(t302)
boxplot(t302)

###(2)

v1<- powerTransform(t302[,1]) ## individual v1
v2<- powerTransform(t302[,2]) ## individual v2
bc<- powerTransform(t302) ## simultaneous about v1 and v2

v1.t302<- bcPower(t302[,1],v1$lambda)
v2.t302<- bcPower(t302[,2],v2$lambda)
bc.t302<- bcPower(t302,bc$lambda)

###individual

par(mfrow=c(1,2))
qqnorm(v1.t302) ## qq graphic for v1
qqnorm(v2.t302) ## qq graphic for v2
par(mfrow=c(1,1))

indit302<- cbind(v1.t302,v2.t302) ## combind v1 and v2
D1<-mahalanobis(indit302,colMeans(indit302),cov(indit302))
qqplot(qchisq(ppoints(25), df = 2), D1,
main = expression("Q-Q plot for" ~~ {chi^2}[nu == 2]))
abline(c(0,1)) ## chi-squre qq plot

f1<-kde2d(indit302[,1],indit302[,2], h = c(width.SJ(indit302[,1]),
width.SJ(indit302[,2],method="dpi")))
persp(f1, phi = 30, theta = 20, d = 5) ## nonparametric

###simultaneous

par(mfrow=c(1,2))
qqnorm(bc.t302[,1]) ## qq graphic
qqnorm(bc.t302[,2]) ## qq graphic
par(mfrow=c(1,1))

D2<-mahalanobis(bc.t302,colMeans(bc.t302),cov(bc.t302))
qqplot(qchisq(ppoints(25), df = 2), D2,
main = expression("Q-Q plot for" ~~ {chi^2}[nu == 2]))
abline(c(0,1)) ## chi-squre qq plot

f2<-kde2d(bc.t302[,1],bc.t302[,2], h = c(width.SJ(bc.t302[,1]),
width.SJ(bc.t302[,2],method="dpi")))
persp(f2, phi = 30, theta = 20, d = 5) ## nonparametric


