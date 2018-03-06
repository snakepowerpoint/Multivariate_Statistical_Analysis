setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####11.7
###(a)

f1 <- function(x){
if (x<=1 & x>=-1)
print(1-abs(x))
else stop("input should be in [-1,1]")
}

f2 <- function(x){
if (x<=1.5 & x>=-0.5)
print((1-abs(x-0.5))/2)
else stop("input should be in [-0.5,1.5]")
}

plot(c(), xlim=c(-2,2), ylim=c(0,2), xlab="x", ylab="f(x)")
curve(f1, xlim=c(-1,1), add=T, col="red", lwd=2, labels="f1")
curve(f2, xlim=c(-0.5,1.5), add=T, col="blue", lwd=2)
legend("topleft", c("f1","f2"), text.col=c(2,4))

###(b)

## R1: f1/f2 >= 1, R2: f1/f2 < 1

###(c)

## R1: f1/f2 >= 0.8/0.2 =4 , R2: f1/f2 < 4


#####7.24
###(a)

t1104 <- read.table("T11-4.dat")
names(t1104)
names(t1104) <- c("X1=CF/TD","X2=NI/TA","X3=CA/CL","X4=CA/NS","Group")

plot(t1104[,1:2], type="n")
text(t1104[,1], t1104[,2], col=c(rep(1,21),rep(2,25)), labels=t1104[,5])
plot(t1104[,c(1,3)], type="n")
text(t1104[,1], t1104[,3], col=c(rep(1,21),rep(2,25)), labels=t1104[,5])
plot(t1104[,c(1,4)], type="n")
text(t1104[,1], t1104[,4], col=c(rep(1,21),rep(2,25)), labels=t1104[,5])
## it seems to exist bivariate normality

###(b)

xbar1 <- apply(t1104[1:21,1:2], 2, mean)
S1 <- cov(t1104[1:21,1:2])
xbar2 <- apply(t1104[21:46,1:2], 2, mean)
S2 <- cov(t1104[21:46,1:2])

library(xtable)
xtable(as.matrix(xbar1))
xtable(as.matrix(S1))
xtable(as.matrix(xbar2))
xtable(as.matrix(S2))

###(c)

k <- log(det(S1)/det(S2))/2 + (xbar1%*%solve(S1)%*%xbar1 - 
xbar2%*%solve(S2)%*%xbar2)/2

## refer to text (11-29)

###(d)

library(MASS)
install.packages("mclust")
library(mclust)

qda.1104 <- qda(t1104[,1:2], t1104[,5], method="mle")
qda.1104
attributes(qda.1104)
pcl <- predict(qda.1104,t1104[,1:2])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER

qda.1104.cv <- qda(t1104[,1:2], t1104[,5], CV=TRUE)
qda.1104.cv$class
table(t1104$Group,qda.1104.cv$class)
mean(t1104$Group!=qda.1104.cv$class) ## E(AER)

###(e)

qda.1104 <- qda(t1104[,1:2], t1104[,5], prior=c(0.05,0.95), method="mle")
qda.1104
pcl <- predict(qda.1104,t1104[,1:2])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER

qda.1104.cv <- qda(t1104[,1:2], t1104[,5], prior=c(0.05,0.95), CV=TRUE)
qda.1104.cv$class
table(t1104$Group,qda.1104.cv$class)
mean(t1104$Group!=qda.1104.cv$class) ## E(AER)
## the prior rates are unreasonable because we get a set of observations
## which p1=21/46, p2=25/46. Moreover, the error rate is high

###(f)

Sp <- ((21-1)*S1 + (25-1)*S2)/(46-2)
(xbar1 - xbar2)%*%solve(Sp)
xtable(Sp)

lda.1104 <- lda(t1104[,1:2], t1104[,5], method="mle")
lda.1104
plot(lda.1104,abbrev=T, col=c(rep(1,21), rep(2,25)))
pcl <- predict(lda.1104, t1104[,1:2])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER
## not a good decision because the error rate is higher that (d)
## that is, quadratic classification rule

###(g)
### (x1,x3)

xbar1 <- apply(t1104[1:21,c(1,3)], 2, mean)
S1 <- cov(t1104[1:21,c(1,3)])
xbar2 <- apply(t1104[21:46,c(1,3)], 2, mean)
S2 <- cov(t1104[21:46,c(1,3)])

xtable(as.matrix(xbar1))
xtable(as.matrix(S1))
xtable(as.matrix(xbar2))
xtable(as.matrix(S2))

k <- log(det(S1)/det(S2))/2 + (xbar1%*%solve(S1)%*%xbar1 - 
xbar2%*%solve(S2)%*%xbar2)/2
## refer to text (11-29)

qda.1104 <- qda(t1104[,c(1,3)], t1104[,5], method="mle")
qda.1104
pcl <- predict(qda.1104,t1104[,c(1,3)])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER

qda.1104.cv <- qda(t1104[,c(1,3)], t1104[,5], CV=TRUE)
qda.1104.cv$class
table(t1104$Group,qda.1104.cv$class)
mean(t1104$Group!=qda.1104.cv$class) ## E(AER)

qda.1104 <- qda(t1104[,c(1,3)], t1104[,5], prior=c(0.05,0.95), method="mle")
qda.1104
pcl <- predict(qda.1104,t1104[,c(1,3)])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER

qda.1104.cv <- qda(t1104[,c(1,3)], t1104[,5], prior=c(0.05,0.95), CV=TRUE)
qda.1104.cv$class
table(t1104$Group,qda.1104.cv$class)
mean(t1104$Group!=qda.1104.cv$class) ## E(AER)
## the prior rates are unreasonable because we get a set of observations
## which p1=21/46, p2=25/46. Moreover, the error rate is high

### (x1,x4)

xbar1 <- apply(t1104[1:21,c(1,4)], 2, mean)
S1 <- cov(t1104[1:21,c(1,4)])
xbar2 <- apply(t1104[21:46,c(1,4)], 2, mean)
S2 <- cov(t1104[21:46,c(1,4)])

xtable(as.matrix(xbar1))
xtable(as.matrix(S1))
xtable(as.matrix(xbar2))
xtable(as.matrix(S2))

k <- log(det(S1)/det(S2))/2 + (xbar1%*%solve(S1)%*%xbar1 - 
xbar2%*%solve(S2)%*%xbar2)/2
## refer to text (11-29)

qda.1104 <- qda(t1104[,c(1,4)], t1104[,5], method="mle")
qda.1104
pcl <- predict(qda.1104,t1104[,c(1,4)])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER

qda.1104.cv <- qda(t1104[,c(1,4)], t1104[,5], CV=TRUE)
qda.1104.cv$class
table(t1104$Group,qda.1104.cv$class)
mean(t1104$Group!=qda.1104.cv$class) ## E(AER)

qda.1104 <- qda(t1104[,c(1,4)], t1104[,5], prior=c(0.05,0.95), method="mle")
qda.1104
pcl <- predict(qda.1104,t1104[,c(1,4)])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER

qda.1104.cv <- qda(t1104[,c(1,4)], t1104[,5], prior=c(0.05,0.95), CV=TRUE)
qda.1104.cv$class
table(t1104$Group,qda.1104.cv$class)
mean(t1104$Group!=qda.1104.cv$class) ## E(AER)
## the prior rates are unreasonable because we get a set of observations
## which p1=21/46, p2=25/46. Moreover, the error rate is high

###(h)

xbar1 <- apply(t1104[1:21,1:4], 2, mean)
S1 <- cov(t1104[1:21,1:4])
xbar2 <- apply(t1104[21:46,1:4], 2, mean)
S2 <- cov(t1104[21:46,1:4])

xtable(as.matrix(xbar1))
xtable(as.matrix(S1))
xtable(as.matrix(xbar2))
xtable(as.matrix(S2))

k <- log(det(S1)/det(S2))/2 + (xbar1%*%solve(S1)%*%xbar1 - 
xbar2%*%solve(S2)%*%xbar2)/2
## refer to text (11-29)

qda.1104 <- qda(t1104[,1:4], t1104[,5], method="mle")
qda.1104
pcl <- predict(qda.1104,t1104[,1:4])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER

qda.1104.cv <- qda(t1104[,1:4], t1104[,5], CV=TRUE)
qda.1104.cv$class
table(t1104$Group,qda.1104.cv$class)
mean(t1104$Group!=qda.1104.cv$class) ## E(AER)

qda.1104 <- qda(t1104[,1:4], t1104[,5], prior=c(0.05,0.95), method="mle")
qda.1104
pcl <- predict(qda.1104,t1104[,1:4])$class
table(pcl)
classError(pcl,t1104[,5]) ## APER

qda.1104.cv <- qda(t1104[,1:4], t1104[,5], prior=c(0.05,0.95), CV=TRUE)
qda.1104.cv$class
table(t1104$Group,qda.1104.cv$class)
mean(t1104$Group!=qda.1104.cv$class) ## E(AER)
## the prior rates are unreasonable because we get a set of observations
## which p1=21/46, p2=25/46. Moreover, the error rate is high


#####11.35
###(a)

t1110 <- read.table("T11-10.dat")
names(t1110)
names(t1110) <- c("Gender","Age","TailLength","BodyLength")

p <- as.numeric(factor(t1110[,2], labels=c(1,2,3)))
plot(t1110[,3:4], col=as.numeric(t1110[,1]), pch=p)
legend("topleft", c("Female Age2","Female Age3","Female Age4","Male Age2", 
"Male Age3","Male Age4"), col=c(1,1,1,2,2,2), pch=c(1,2,3,1,2,3), 
text.col=c(1,1,1,2,2,2))
## we can separate these snakes by gender and age

###(b)

lda.1110 <- lda(t1110[,3:4], t1110[,1], CV=TRUE)
table(lda.1110$class)
classError(lda.1110$class,t1110[,1]) ## E(AER)

###(c)

lda.1110 <- lda(t1110[,3:4], t1110[,2], CV=TRUE)
table(lda.1110$class)
classError(lda.1110$class,t1110[,2]) ## E(AER)

###(d)

lda.1110 <- lda(Age ~ BodyLength, t1110, CV=TRUE, method="mle")
table(lda.1110$class)
classError(lda.1110$class,t1110[,2])
## cannot use only one variable to separate observations


