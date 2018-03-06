setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####6.24
library(xtable)

t613 <- read.table("T6-13.dat")
library(car)

t613$V5 <- as.factor(t613$V5)
fit.lm <- lm(cbind(V1,V2,V3,V4)~ V5 , data = t613)
summary(fit.lm)

bone.manova<- Manova(fit.lm)
summary(bone.manova, test="Wilks") 
## use Wilks test
## lambda=0.8301, F=2.049, P(F(8,168)>2.049)=0.044, so reject


trteff.diff<-function (data, level=0.95) {
# Bonferroni Based Simultaneous CI for Treatments Difference for one-way MANOVA
# data: Nx(p+1) numeric matrix or dataframe of data, N is the total sample size, p  is number of variables,
# the last column is the factors
# level: confidence level of interval,default=0.95

    Y<-data[,-ncol(data)]
    X<-data[,ncol(data)]
    p <- ncol(Y) 
    g <- length(levels(as.factor(X)))
    X <- as.numeric(X)
    data <- as.matrix(cbind(Y,X))
     
    N <- length(X)
    meanvec <- matrix(apply(Y, 2, mean), ncol = 1)
    n <- matrix(rep(0, g), ncol = 1)
    trtmean <- matrix(rep(0, p * g), ncol = g)
    trt.effect <- matrix(rep(0, p * g), ncol = g)
    trt.cov <- matrix(rep(0, (p * p) * g), nrow = p)
    W <- matrix(rep(0, (p * p)), nrow = p)
 
    df2 <- 0
    for (k in 1:g) {
       n[k] <- length(subset(X, X == k))
       trtmean[, k] <- as.matrix(colMeans(subset(Y, X == k)))
       trt.effect[, k] <- trtmean[, k] - meanvec
      for (i in 1:p) {
        for (j in 1:p) {
          trt.cov[j, i + (k - 1) * p] <- cov(subset(Y[,i], X == k), subset(Y[, j], X == k))
         }
       }
    W = W + (n[k] - 1) * trt.cov[, (1 + (k - 1) * p):(k *p)]
   }
    
    prob <- 1 - ((1 - level)/(p * g * (g - 1)))
    t <- qt(prob, (N - g))
    d <- (g * (g - 1)/2)
    tau <- matrix(rep(0, (p * d)), ncol = 1)
    lower <- matrix(rep(0, (p * d)), ncol = 1)
    upper <- matrix(rep(0, (p * d)), ncol = 1)
    C<-matrix(0,nrow=d,g)
    kk<-1
    for(k in 1:(g-1))
     for(l in (k+1):g){
       C[kk,k]<-1
       C[kk,l]<--1
       kk<-kk+1
       }
    for (i in 1:p) {
       for (k in 1:(g - 1)) {
         for (l in (k + 1):g) {      
           tau[(1 + (i - 1) * d):(i * d), ] <- C%*%(as.matrix(trt.effect[i, ]))
           lower[(1 + (i - 1) * d):(i * d), ] <- tau[(1 +(i - 1) * d):(i * d), ] - t * sqrt((W[i, i]/(N -
               g)) * ((1/n[k]) + (1/n[l])))
           upper[(1 + (i - 1) * d):(i * d), ] <- tau[(1 +(i - 1) * d):(i * d), ] + t * sqrt((W[i, i]/(N -
               g)) * ((1/n[k]) + (1/n[l])))
         }
      }
     }
cat(trt.effect,"\n")
cat("\n Bonferroni Based Simultaneous CI for Treatments Difference \n")
SCI = data.frame(Trt = C, Estimate = round(tau,3), LowerCI = round(lower, 3), UpperCI = round(upper,3))
labs<-NA
for(i in 1:p)
 for(k in 1:(g-1))
   for(l in (k+1):g)
     labs<-c(labs,paste("tau",k,"-",l,i,sep=""))
rownames( SCI)<-labs[-1]
print(SCI)
}


trteff.diff(t613)
xtable(trteff.diff(t613))
## all intervals include 0, so we chooese the most possible one
## tau1-31
## tau1-33
## tau2-33

plot(t613)

chiqqplot<-function(x){
	if(!is.matrix(x))
	x<-as.matrix(x)
	n<-nrow(x)
	p<-ncol(x)
	D2<-mahalanobis(x,colMeans(x),cov(x))
	qqplot(qchisq(ppoints(n), df = p), D2,xlab="Theoretical Q of Chisquare", ylab="Mah distance",
       main = expression("Q-Q plot for" ~~ {chi^2}[p]),pch=19)
abline(c(0,1))
}

chiqqplot(t613[,1:4]) ## normal
cov(t613[1:30,1:4])
cov(t613[31:60,1:4])
cov(t613[61:90,1:4]) ## the var-cov matrixs are similar

#####6.25

t1107 <- read.table("T11-7.dat")
factor <- factor(t1107$V6,levels=c("Wilhelm","SubMuli","Upper"),labels=c('1','2','3'))
t1107 <- t1107[,-6]
t1107 <- cbind(t1107,factor)
t1107$factor <- as.factor(t1107$factor)

fit.lm1 <- lm(cbind(V1,V2,V3,V4,V5)~factor, data=t1107)
summary(fit.lm1)

oil.manova <- Manova(fit.lm1)
summary(oil.manova, test="Wilks") 
qf(0.95,10,98)
## use Wilks test
## lambda=0.1159, F=18.9848, P(F(10,98)>18.9848)=0.00..., so reject

trteff.diff(t1107)
xtable(trteff.diff(t1107))

plot(t1107) ## var3 is abnormal
chiqqplot(t1107[,1:5]) ## not normal enough
cov(t1107[1:7,1:5])
cov(t1107[8:18,1:5])
cov(t1107[19:56,1:5]) ## var-cov matrixs are different

## example 11.14
v2.t1107 <- sqrt(t1107[,2])
v3.t1107 <- sqrt(t1107[,3])
v4.t1107 <- t1107[,4]^(-1)
mt1107 <- cbind(t1107$V1,v2.t1107,v3.t1107,v4.t1107,t1107$V5,t1107$factor)
mt1107 <- data.frame(mt1107)
colnames(mt1107) <- c("V1","V2","V3","V4","V5","factor")
mt1107$factor <- as.factor(mt1107$factor)

fit.lm2 <- lm(cbind(V1,V2,V3,V4,V5)~factor, data=mt1107)
summary(fit.lm2)

oil.manova1 <- Manova(fit.lm2)
summary(oil.manova1, test="Wilks")
## use Wilks test
## lambda=0.1198, F=18.5171, P(F(10,98)>18.9848)=0.00..., so reject

chiqqplot(mt1107[,1:5])

#####6.27
###(a)

t614 <- read.table("T6-14.dat")
mt614 <- data.frame(t614[1:30,])
ft614 <- data.frame(t614[31:60,])

mbar <- colMeans(mt614[,1:4])
fbar <- colMeans(ft614[,1:4])
Sp <- ((30-1)*cov(mt614[,1:4])+(30-1)*cov(ft614[,1:4]))/(30+30-2)

plot(fbar, pch=20, col='red')
lines(fbar, col='red')
points(mbar, pch=20, col='blue')
lines(mbar, col='blue')

###(b)

C <- matrix(0,3,4)
for (i in 1:3){
 C[i,i:(i+1)]=c(-1,1)
}

t(C%*%(mbar-fbar))%*%solve(C%*%Sp%*%t(C))%*%C%*%(mbar-fbar)*(1/30+1/30)^(-1)
qf(0.95,4-1,30+30-4)*(30+30-2)*(4-1)/(30+30-4) ## can not reject

I <- matrix(1,4,1)

t(I)%*%(mbar-fbar)%*%solve(t(I)%*%Sp%*%I)%*%t(I)%*%(mbar-fbar)*(1/30+1/30)^(-1)
qf(0.95,1,30+30-2) ## can not reject

xbar <- (30*mbar+30*fbar)/(30+30)
S <- cov(t614[,1:4])

t(C%*%xbar)%*%solve(C%*%S%*%t(C))%*%C%*%xbar
qf(0.95,4-1,30+30-4+1)*(30+30-1)*(4-1)/(30+30-4+1) ## can not reject

#####6.31
###(a)

t617 <- read.table("T6-17.dat", col.names=c("location", "variety", 
"yield", "weight", "seedsize"))

t617$location <- as.factor(t617$location)
t617$variety <- as.factor(t617$variety)

fit.lm3 <- lm(cbind(yield, weight, seedsize) ~ 
location * variety, data = t617)

peanuts.manova <- Manova(fit.lm3)
summary(peanuts.manova, test="Wilks") 
## check location, variety and interaction
## location: P(F(3,4)>11.1843)=0.0205
## variety: P(F(6,8)>10.6191=0.0019
## location*variety: P(F(6,8)>3.5582=0.0508)
## with 95% confident level, we reject location and variety

###(b)

attributes(peanuts.manova)
peanuts.manova$SSPE

means.matrix <- matrix(0,12,3)
for (i in seq(1,11,2)) {
 means.matrix[i,] = colMeans(t617[i:(i+1),3:5])
 means.matrix[i+1,] = colMeans(t617[i:(i+1),3:5])
}

error <- cbind(t617[,1:2],t617[,3:5]-means.matrix)

par(mfrow=c(1,3))
qqnorm(error[,3])
qqline(error[,3]) ## check normality of yield

qqnorm(error[,4])
qqline(error[,4]) ## check normality of weight

qqnorm(error[,5])
qqline(error[,5]) ## check normality of seedsize
par(mfrow=c(1,1))

plot(error[,3:5]) ## check simultaneous normality
chiqqplot(error[,3:5]) ## the residual vectors have outliers

###(c)

fit.lm4 <- lm(yield~ location * variety, data = t617)
fit.lm5 <- lm(weight~ location * variety, data = t617)
fit.lm6 <- lm(seedsize~ location * variety, data = t617)

anova(fit.lm4) ## significant in interaction
anova(fit.lm5) ## significant in interaction
anova(fit.lm6) ## not significant in interaction (seedsize)
xtable(anova(fit.lm4))
xtable(anova(fit.lm5))
xtable(anova(fit.lm6))

###(d)

variety <- cbind(t617[c(3,4,7,8,11,12),3:5],t617[c(3,4,7,8,11,12),2])
colnames(variety) <- c("yield","weight","seedsize","factor")

trteff.diff(variety)
## the sample size is too small to do test, but we can infer that
## variety8 is better than the others in every characteristic
xtable(trteff.diff(variety))

#####nonsense

W <- cov(t617[3:4,3:5]) + cov(t617[7:8,3:5]) + cov(t617[11:12,3:5])

means.matrix[3,3]-means.matrix[11,3]+
sqrt(var(t617[3:4,5])*(1/2+1/2))*qt(1-0.05/(1*3*(3-1)),6-3)

means.matrix[3,3]-means.matrix[11,3]-
sqrt(var(t617[3:4,5])*(1/2+1/2))*qt(1-0.05/(1*3*(3-1)),6-3)


#####
install.packages("quantreg")
library(quantreg)
library(xtable)

xtable(cov(t613[1:30,1:4]))
xtable(cov(t613[31:60,1:4]))
xtable(cov(t613[61:90,1:4]))





