#----------------------------------------------------------------------------
# FLDA 
#----------------------------------------------------------------------------
library(RColorBrewer)
showMatrix <- function(x, ...)
  image(t(x[nrow(x):1,]), xaxt='none', yaxt='none', 
  col=rev(colorRampPalette(brewer.pal(9,'Greys'))(100)), ...)

faces <- read.table("faces.txt")
x <- as.matrix(faces)
rownames(x) <- paste("P", rep(1:40, each=10), "-", rep(1:10,40), sep="")
john <- t(faces)

## Here we take persons as variables. All we want to do is find out 
## the variations which are made by those specific ones.

library(psych)
cortest.bartlett(john)
fa.parallel(john) ## the number of factors=34, components=27

sum(eigen(cor(john))$values[1:15])/sum(eigen(cor(john))$values)
## 78% variation

fa.faces <- factanal(john, factors=15, scores="regression", 
rotation="varimax", method="mle")

fa.faces$loadings 
## factor 1: 231~238 is small. 24th 
## factor 1: 330~339, 104~109. 34th, 11th

load <- fa.faces$loadings
as.matrix(as.numeric(load[,1]>0.7)) ## 25th, 23th, 10th, 9th
as.matrix(as.numeric(load[,2]>0.7)) ## 17th, 5th,
as.matrix(as.numeric(load[,3]>0.7)) ## 11th

plot(fa.faces$scores)
plot(t(fa.faces$scores))

d.faces <- dist(faces)
ha <- hclust(d.faces, method="complete")
windows()
plot(ha, hang=-1)

#----------------------------------------------------------------------------
# Single, Complete, Average toward 5th and 40th persons
#----------------------------------------------------------------------------
xs <- x[c(41:50,391:400),]
xd <- dist(xs)
par(mfrow=c(4,5), mar=rep(0,4))
for(i in 1:20)
showMatrix(matrix(xs[i,],64,64))

library(TeachingDemos)
windows()
ha <- hclust(xd, method="average")
plot(ha, hang=-1)
for(i in 1:nrow(xs))
subplot(showMatrix(matrix(xs[ha$order[i],],64,64)), i, 70, size=c(0.5,0.5))

windows()
hc <- hclust(xd, method="complete")
plot(hc, hang=-1)
for(i in 1:nrow(xs))
subplot(showMatrix(matrix(xs[hc$order[i],],64,64)), i, 70, size=c(0.5,0.5))

windows()
hs <- hclust(xd, method="single")
plot(hs, hang=-1)
for(i in 1:nrow(xs))
subplot(showMatrix(matrix(xs[hs$order[i],],64,64)), i, 70, size=c(0.5,0.5))

## on the other hand, 5th and 6th can separate

xss <- x[c(41:50,51:60),]
xdd <- dist(xss)

haa <- hclust(xdd, method="average")
plot(haa, hang=-1)
for(i in 1:nrow(xss))
subplot(showMatrix(matrix(xss[haa$order[i],],64,64)), i, 70, size=c(0.5,0.5))


#----------------------------------------------------------------------------
# cluster 10 persons
#----------------------------------------------------------------------------
xss <- x[c(seq(10,400,by=10)),]
par(mfrow=c(4,10), mar=rep(0,4))
for(i in 1:40)
showMatrix(matrix(xss[i,],64,64))

xdd <- dist(xss)
hh <- hclust(xdd, method="average")
plot(hh, hang=-1)

#----------------------------------------------------------------------------
# k-nn
#----------------------------------------------------------------------------
library(class)

knn.faces <- knn(train=ft[,-1], test=fv[,-1], cl=ft[,1], k=40)
knn.faces
table(knn.faces)
classError(knn.faces, fv[,1])

knn.faces.1 <- knn(train.faces[,-1], test.faces[,-1], cl=train.faces[,1], 
k=40)
knn.faces.1
table(knn.faces.1)
classError(knn.faces.1, test.faces[,1])

#----------------------------------------------------------------------------
# qda
#----------------------------------------------------------------------------
qda.faces <- qda(ft[,-1], ft[,1], method="mle") ## useless

#----------------------------------------------------------------------------
# SVM
#----------------------------------------------------------------------------
library(kernlab)
library(xtable)

clist = 2^seq(-4,4)
err <- matrix(0,9,2)
for (i in seq(length(clist))){
  svm <- ksvm(X1~., data=ft, type="C-svc", kernel="polydot", 
          kpar=list(degree=2), C=clist[i], cross=5)
  pre.1 <- predict(svm, ft[,-1])
  pre.2 <- predict(svm, fv[,-1])
  err[i,1] <- classError(pre.1, ft[,1])$errorRate
  err[i,2] <- classError(pre.2, fv[,1])$errorRate
}
colnames(err) = c("train", "test")
print(cbind(clist,err))
## for any C, the error rate is 0.1 as degree=2
xtable(print(cbind(clist,err)))

err.1 <- matrix(0,9,2)
for (i in seq(length(clist))){
  svm <- ksvm(X1~., data=ft, type="C-svc", kernel="polydot", 
          kpar=list(degree=1), C=clist[i], cross=5)
  pre.1 <- predict(svm, ft[,-1])
  pre.2 <- predict(svm, fv[,-1])
  err.1[i,1] <- classError(pre.1, ft[,1])$errorRate
  err.1[i,2] <- classError(pre.2, fv[,1])$errorRate
}
colnames(err.1) = c("train", "test")
print(cbind(clist,err.1))
## for any C, the error rate is 0 as degree=1
xtable(print(cbind(clist,err.1)))

svm.faces <- ksvm(X1~., data=ft, type="C-svc", kernel="polydot", 
kpar=list(degree=1), C=1, cross=5)
presvm <- predict(svm.faces, ft[,-1])
classError(presvm, ft[,1]) ## train set
classError(predict(svm.faces, fv[,-1]), fv[,1])
predict(svm.faces, fv[,-1])

#----------------------------------------------------------------------------

svm.faces.1 <- ksvm(X1~., data=ft, type="C-svc", kernel="laplacedot", 
kpar=list(sigma=0.05), C=1, cross=5)
presvm.1 <- predict(svm.faces.1, ft[,-1])
classError(presvm.1, ft[,1]) ## train set
classError(predict(svm.faces.1, fv[,-1]), fv[,1])





#############################################################################
#           Demonstration
#############################################################################
## load data

library(RColorBrewer)
showMatrix <- function(x, ...)
  image(t(x[nrow(x):1,]), xaxt='none', yaxt='none', 
  col=rev(colorRampPalette(brewer.pal(9,'Greys'))(100)), ...)

faces <- read.table("faces.txt")
dim(faces)
x <- as.matrix(faces)
par(mfrow=c(20,20), mar=rep(0,4))
for(i in 1:400)
showMatrix(matrix(x[i,],64,64))
rownames(x) <- paste("P", rep(1:40, each=10), "-", rep(1:10,40), sep="")

## take 39th for example

xs <- x[381:390,]
xd <- dist(xs)
par(mfrow=c(2,5), mar=rep(0,4))
for(i in 1:10)
showMatrix(matrix(xs[i,],64,64))

## use Hierarchical clustering and the method by Bien and Tibshirani(2011)
## to cluster the data

library(TeachingDemos)

par(mfrow=c(2,1))
hh <- hclust(xd, method="average")
plot(hh, hang=-1)
for(i in 1:nrow(xs))
subplot(showMatrix(matrix(xs[hh$order[i],],64,64)), i, 70, size=c(0.5,0.5))
## perform minimax linkage clustering by Bien, J., and Tibshirani, R.(2011)

install.packages("protoclust")
library(protoclust)
hm <- protoclust(xd)
plot(hm, hagn=-1)
for(i in 1:nrow(xs))
subplot(showMatrix(matrix(xs[hm$order[i],],64,64)), i, 70, size=c(0.5,0.5))

## compute the coordinates of nodes
absi <- function(hc, level=length(hc$height), init=TRUE){
  if(init){
    .cout <<- 0
    .topAbsis <<- NULL
    .heights <<- NULL
  }
  if(level<0){
     .count <<- .count + 1
     return(.count)
  }
  node <- hc$merge[level,]
  le <- absi(hc, node[1], init=FALSE)
  ri <- absi(hc, node[2], init=FALSE)
  mid <- (le+ri)/2
  .topAbsis <<- c(.topAbsis, mid)
  .heights <<- c(.heights, hc$height[level])
  invisible(mid)
}

absi(hm)
ord <- order(.heights)
yo <- .heights[ord]
xo <- .topAbsis[ord]
for(i in 1:(nrow(xs)-1))
subplot(showMatrix(matrix(xs[hm$proto[i],64,64)),xo[i],yo[i],size=c(0.5,0.5))


#----------------------------------------------------------------------------
#         Discrimination
#----------------------------------------------------------------------------
id <- rep(1:40, each=10)
faces.data.frame <- data.frame(cbind(id=id, faces))
dim(faces.data.frame)
## pick out the last photo for each person as the testing photo

testid <- seq(10, 400, by=10)
train.faces <- faces.data.frame[-testid,]
test.faces <- faces.data.frame[testid,]

## calculate the directions of principal components, determine the number 
## of principal components, and calculate thier scores

xc <- scale(train.faces[,-1], scale=FALSE)
A <- t(xc)/sqrt(360-1)
dim(A) 
## 4096x360
## thus, the covariance matrix is A%*%t(A)
A.egn <- eigen(t(A)%*%A)
## 360x360
pc <- A%*%A.egn$vectors
pc <- apply(pc,2,function(i) i/sqrt(sum(i*i)))
## normalize the pc
n <- 80
sum(A.egn$value[1:n])/sum(A.egn$value)
## 92%, thus we use the first 80 pcs
pcs <- pc[,1:n]
## pc scores for training data
yt <- xc%*%pcs

plot(yt, type='n', xlab="Component 1", ylab="Component 2")
text(yt[,1], yt[,2], labels=as.matrix(id)[-testid,], col=
as.matrix(id)[-testid,])

plot(yt[,c(1,3)], type='n', xlab="Component 1", ylab="Component 3")
text(yt[,1], yt[,3], labels=as.matrix(id)[-testid,], col=
as.matrix(id)[-testid,])

plot(yt[,c(2,3)], type='n', xlab="Component 2", ylab="Component 3")
text(yt[,2], yt[,3], labels=as.matrix(id)[-testid,], col=
as.matrix(id)[-testid,])

## use discrimination function under PC space, and we find an excellent
## result

ft <- data.frame(cbind(id[-testid], yt))
head(ft)
library(MASS)
flda <- lda(X1~.,ft)
fl <- predict(flda, ft)$class
diag(table(fl, id[-testid])) ## the result is perfect!!
table(fl)
## cannot know whether 1 classified into 2, and 2 classified into 1

fl
id[-testid]
table(fl, id[-testid])
library(mclust)
classError(fl, id[-testid])

## calculate the scores of the test groups under PC space,
## use the discrimination function above to classify the test groups.
## we will get a nearly perfect result

##pc scores for testing data
xv <- as.matrix(test.faces[,-1])
xvs <- scale(xv, scale=FALSE)
yv <- xvs%*%pcs
fv <- data.frame(cbind(id[testid], yv))
pl <- predict(flda, fv)$class
diag(table(pl, id[testid])) ## number 5 misses
table(pl)

pl
id[testid]
classError(pl, id[testid])
## number 5 and 40 are similar


#----------------------------------------------------------------------------
#         Appendix
#----------------------------------------------------------------------------
kmo <- function(x)
{
x <- subset(x, complete.cases(x)) # Omit missing values
r <- cor(x) # Correlation matrix
r2 <- r^2 # Squared correlation coefficients
i <- solve(r) # Inverse matrix of correlation matrix
d <- diag(i) # Diagonal elements of inverse matrix
p2 <- (-i/sqrt(outer(d, d)))^2 # Squared partial correlation coefficients
diag(r2) <- diag(p2) <- 0 # Delete diagonal elements
KMO <- sum(r2)/(sum(r2)+sum(p2))
MSA <- colSums(r2)/(colSums(r2)+colSums(p2)) #individual measures of sampling adequacy for each item
return(list(KMO=KMO, MSA=MSA))
}




