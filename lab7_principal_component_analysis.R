setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####exercise1
#####8.26

###using S
###(a)

t406 <- read.table("T4-6.dat")
names(t406)
names(t406) <- c("indep","supp","benev","conform","leader","gen","posit")

pc.406 <- princomp(t406[,1:5])
summary(pc.406, loadings=T)
plot(pc.406, type="lines") ## choose 1~4

###(b)

library(xtable)
load.406 <- matrix(loadings(pc.406),5,5)
load.406
xtable(load.406)

## component1 represents difference between indep and benev & conform
## component2 represents difference between supp and comform & leader
## component3 represents difference between leader and indep & conform
## component4 represents difference between supp and benev

###(c)

plot(pc.406$scores[,1:2], type="n")
text(pc.406$scores[,1], pc.406$scores[,2], labels=t406$gen, cex=0.7, lwd=2,
col=c(t406$gen))
## there are some outliers

text(pc.406$scores[,1], pc.406$scores[,2], labels=1:130, cex=0.7, lwd=2)

###(d)

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

evalue.CI(t406[,1:5],1)
eigen(cov(t406[,1:5])) ## notice that the eigenvalue is based on S
130*pc.406$sdev^2/129 ## unbiased eigenvalues

pc.406.a <- prcomp(t406[,1:5])
summary(pc.406.a)
summary(pc.406)
pc.406.a$sdev^2 ## same as above

###using R
###(a)

pc.406.r <- princomp(t406[,1:5], cor=T)
summary(pc.406.r, loadings=T)
summary(pc.406, loadings=T) ## a little bit different
plot(pc.406.r, type="lines") ## also choose 1~4

###(b)

load.406.r <- matrix(loadings(pc.406.r),5,5)
load.406.r
xtable(load.406.r)
## similar with results from S

###(c)

plot(pc.406.r$scores[,1:2], type="n")
text(pc.406.r$scores[,1], pc.406.r$scores[,2], labels=t406$gen, cex=0.7,
lwd=2, col=c(t406$gen))
## there are some outliers

text(pc.406.r$scores[,1], pc.406.r$scores[,2], labels=1:130, cex=0.7, lwd=2)

#####exercise2

install.packages("plsgenomics")
library(plsgenomics)
data(leukemia)
names(leukemia)
str(leukemia)
attributes(leukemia)

leukemia$Y
head(leukemia$X)
leukemia$gene.names[,3]
leukemia$gene.names[3051,] ## y is description of gene
dim(leukemia$gene.names)
dim(leukemia$X) ## x is data

pc.leu <- princomp(leukemia$X)
pc.leu <- prcomp(leukemia$X)
summary(pc.leu)
attributes(pc.leu)
str(pc.leu)

plot(1:38, pc.leu$sdev^2, type="b", xlab="Components", ylab="eigenvalue")
plot(pc.leu, type="lines") ## choose the first 10 components

dim(pc.leu$rotation) ## eigen vactors
head(pc.leu$rotation)
dim(pc.leu$x) ## x represents values of the first 38 components

eigen(cov(leukemia$X))$values 
## some are very small
## so we can use first 38 components

pc.leu$x[,1:10]
plot(pc.leu$x[,1:2], type="n")
text(pc.leu$x[,1], pc.leu$x[,2], labels=c(rep(1,27),rep(2,11)),
col=c(rep(1,27),rep(2,11)))

par(mfrow=c(3,3))

for (i in 2:10){
plot(pc.leu$x[,c(1,i)], type="n")
text(pc.leu$x[,1], pc.leu$x[,i], labels=c(rep(1,27),rep(2,11)),
col=c(rep(1,27),rep(2,11)))
i=i+1
}

plot(pc.leu$x[,1], type="n")
text(pc.leu$x[,1],labels=c(rep(1,27),rep(2,11)),
col=c(rep(1,27),rep(2,11)))

pc.leu$rotation[,1]

#----------------------------------------------------------------
#####exercise1

# Pricipal Components Analysis
# entering raw data and extracting PCs
# from the correlation matrix
mydata<-iris[,1:4]
fit1 <- princomp(mydata, cor=TRUE) #variances are computed with the divisor N
fit2 <- prcomp(mydata,cor=TRUE) #variances are computed with the usual divisor N-1
summary(fit1) # print variance accounted for
loadings(fit1) # pc loadings
plot(fit1,type="lines") # scree plot
fit1$scores # the principal components
biplot(fit1)

#####exercise2

install.packages("plsgenomics")
library(plsgenomics)
data(leukemia)
?leukemia
names(leukemia)
