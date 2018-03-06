setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####9.19
###(a)

t912 <- read.table("T9-12.dat")
names(t912)
names(t912) <- c("SaleGrow","SaleInt","NewSale","Create","Machine",
"Abstract","Math")

t912.s <- scale(t912, center=T, scale=T)

library(xtable)

## m=2, mle
t912.fa2 <- factanal(t912.s, factors=2, scores="regression",
rotation="none", method="mle")
t912.fa2

john <- xtable(t(matrix(t912.fa2$loadings,7,2)))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

## m=3, mle
t912.fa3 <- factanal(t912.s, factors=3, scores="regression",
rotation="none", method="mle")
t912.fa3

john <- xtable(t(matrix(t912.fa3$loadings,7,3)))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

###(b)

## m=2, mle, with rotation
t912.fa2 <- factanal(t912.s, factors=2, scores="regression",
rotation="varimax", method="mle")
t912.fa2
## factor1 may be the ability of math
## factor2 may be the ability of creation

john <- xtable(t(matrix(t912.fa2$loadings,7,2)))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

## m=3, mle, with rotation
t912.fa3 <- factanal(t912.s, factors=3, scores="regression",
rotation="varimax", method="mle")
t912.fa3
## factor1 may be the ability of math
## factor2 may be the ability of creation
## factor3 may be the ability of abstract

john <- xtable(t(matrix(t912.fa3$loadings,7,3)))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

###(c)

t912.fa2
attributes(t912.fa2)
1-t912.fa2$unique
t912.fa2$unique
es2 <- t912.fa2$loadings%*%t(t912.fa2$loadings) + diag(t912.fa2$unique)
cor(t912) - es2

john <- xtable(matrix(1-t912.fa2$unique,1,7))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

john <- xtable(diag(t912.fa2$unique))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

john <- xtable(matrix(es2,7,7))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

t912.fa3
1-t912.fa3$unique
t912.fa3$unique
es3 <- t912.fa3$loadings%*%t(t912.fa3$loadings) + diag(t912.fa3$unique)
cor(t912) - es3
## m=3 is more closer to the matrix R

john <- xtable(matrix(1-t912.fa3$unique,1,7))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

john <- xtable(diag(t912.fa3$unique))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

john <- xtable(matrix(es3,7,7))
print(john, floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

###(d)

## m=2
(50-1-(2*7+4*2+5)/6)*log(det(es2)/det(cor(t912)))
qchisq(0.99, ((7-2)^2-7-2)/2) ## significant

## m=3
(50-1-(2*7+4*3+5)/6)*log(det(es3)/det(cor(t912)))
qchisq(0.99, ((7-3)^2-7-3)/2) ## significant
## choose model with m=3

###(e)

x <- matrix(c(110, 98, 105, 15, 18, 20, 35),1,7)
x.s <- scale(x, center=colMeans(t912), scale=sqrt(diag(cov(t912))))
loading.912 <- matrix(t912.fa3$loadings,7,3)
uniq.912 <- diag(t912.fa3$uniq)

## WLS
wls <- solve(t(loading.912)%*%solve(uniq.912)%*%loading.912)%*%
t(loading.912)%*%solve(uniq.912)%*%t(x.s)
wls

print(xtable(wls), floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

## regression

reg <- t(loading.912)%*%solve(cor(t912))%*%t(x.s)
reg

print(xtable(reg), floating=FALSE, tabular.environment="array",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

#####9.28

t109 <- read.table("T1-9.dat")
names(t109)
names(t109) <- c("Country", "100", "200", "400", "800.m", "1500.m",
"3000.m", "Marathon.m")

install.packages("psych")
library(psych)

fa.parallel(t109[,2:8]) ## numbers of factor are 2

cov.109 <- matrix(cov(t109[,2:8]),7,7)

library(GPArotation)

## use matrix S, m=2, mle
t109.fa.s <- fa(t109[,2:8], covar=T, nfactors=2, scores="regression",
rotate="varimax", fm="ml")
t109.fa.s
plot(t109.fa.s$scores, type="n")
text(t109.fa.s$scores[,1], t109.fa.s$scores[,2], labels=1:nrow(t109))
## observation 11,40 are outliers

## use matrix R, m=2, mle
t109.fa.r <- fa(t109[,2:8], nfactors=2, scores="regression",
rotate="varimax", fm="ml")
t109.fa.r
plot(t109.fa.r$scores, type="n")
text(t109.fa.r$scores[,1], t109.fa.r$scores[,2], labels=1:nrow(t109))
## observation 11,31,46 are outliers

cov(t109[,2:8])
## the results built on S and R respectively are different because the 
## cov is dominated by variable Marathon which is recorded in minutes



