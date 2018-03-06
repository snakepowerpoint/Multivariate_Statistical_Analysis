setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#-------------------------------------------------------------------
digits<-read.table("digits.txt")
digits<-as.matrix(digits)

library(RColorBrewer)

showMatrix <- function(x, ...) image(t(x[nrow(x):1,]), xaxt = 'none', 
yaxt = 'none', col = rev(colorRampPalette(brewer.pal(9, 'GnBu'))(100)), ...)

par(mfrow=c(1,3),mar=rep(0,4))
showMatrix(matrix(digits[1,],16,16))
showMatrix(matrix(digits[1101,],16,16))
showMatrix(matrix(digits[2201,],16,16))

#let¡¦s start with lda classifier.
library(MASS)
digits<-data.frame(cl=rep(1:3,each=1100),digits)
#take 75% samples as training set, you can change this
idx<-c(1:825,1100+1:825,2200+1:825)
#you should take the training set randomly.
train<-digits[idx,]
## record the running time
ptm <- proc.time()
z<-lda(cl~.,data=train)
proc.time() - ptm
plot(z, col=train[,1])

#Now it¡¦s your turn to complete...
#-------------------------------------------------------------------

###(1)

install.packages("mclust")
library(mclust)

tr1 <- sample(1:1100, 825)
tr2 <- sample(1101:2200, 825)
tr3 <- sample(2201:3300, 825)
train <- rbind(digits[tr1,], digits[tr2,], digits[tr3,])
test <- digits[-c(tr1,tr2,tr3),]

###(2)&(3)
###flda

ptm1 <- proc.time()
flda.digits <- lda(cl~., data=train)
proc.time() - ptm1 ## (2.08  0.03  3.16)
attributes(flda.digits)
plot(flda.digits, col=train[,1])
pc1 <- predict(flda.digits, test[,-1])$class
table(pc1)
classError(pc1,test[,1]) ## 0.025

install.packages("kernlab")
library("kernlab")

###svm

# Train linear SVM 
clist <- 2^seq(-6,6)
errlin <- numeric(length(clist))

for (i in seq(length(clist))) {
  svp <- ksvm(cl~.,data=train,type="C-svc",kernel='vanilladot',C=clist[i],scaled=c(),cross=3)
  errlin[i] <- cross(svp)
}

# Plot the CV error as a function of C
plot(clist,errlin,type='l',log="x",ylim=c(0,1),xlab="C",ylab="Error rate")
grid()

which.min(errlin) ## 7th
min(errlin)
clist[7] ## 1

ptm2 <- proc.time()
svm.digits <- ksvm(cl~., data=train, type="C-svc", kernel="vanilladot", 
C=1, scaled=c(),cross=3)
proc.time() - ptm2 ## (2.14  0.05  3.54)

predict(svm.digits, test[,-1], type="decision")
predict(svm.digits, test[,-1], type="vote")
pc2 <- predict(svm.digits, test[,-1])
table(pc2)
classError(pc2,test[,1]) ## 0.025

###kNN

library(class)

cl <- train[,1]
ptm3 <- proc.time()
knn.digits <- knn(train[,-1], test[,-1], cl, k=5, prob=T, use.all=T)
proc.time() - ptm3 ## (5.99  0.10  7.72)
attributes(knn.digits)

table(knn.digits)
classError(knn.digits, test[,1]) ## 0.0048

###(4)

prin.digits <- princomp(digits[,-1])
summary(prin.digits, loadings=T)
summary(prin.digits)
screeplot(prin.digits)
plot(1:length(prin.digits$sdev), prin.digits$sdev^2, type="b",
main="Scree Plot", xlab="Number of Components", ylab="Eigenvalue")

install.packages("psych")
library(psych)

fa.parallel(digits[,-1])
## suggested numbers of factor are 40, numbers of component are 36

digits.1 <- prin.digits$scores[,1:36]
plot(digits.1[,c(1,2)], col=digits[,1])

digits.1 <- data.frame(cl=rep(1:3,each=1100),digits.1)
train <- rbind(digits.1[tr1,], digits.1[tr2,], digits.1[tr3,])
test <- digits.1[-c(tr1,tr2,tr3),]
dim(train)
dim(test)

ptm1 <- proc.time()
flda.digits <- lda(cl~., data=digits.1)
proc.time() - ptm1 ## (0.27  0.01  1.3)
plot(flda.digits, col=train[,1])
pc1 <- predict(flda.digits, test[,-1])$class
table(pc1)
classError(pc1,test[,1]) ## 0.0206

###svm

# Train linear SVM 
clist <- 2^seq(-6,6)
errlin <- numeric(length(clist))

for (i in seq(length(clist))) {
  svp <- ksvm(cl~.,data=train,type="C-svc",kernel='vanilladot',C=clist[i],scaled=c(),cross=3)
  errlin[i] <- cross(svp)
}

# Plot the CV error as a function of C
plot(clist,errlin,type='l',log="x",ylim=c(0,1),xlab="C",ylab="Error rate")
grid()

which.min(errlin) ## 7th
min(errlin)
clist[7] ## 1

ptm2 <- proc.time()
svm.digits <- ksvm(cl~., data=train, type="C-svc", kernel="vanilladot", 
C=1, scaled=c(),cross=3)
proc.time() - ptm2 ## (117.11  0.03  118.72)

predict(svm.digits, test[,-1], type="decision")
predict(svm.digits, test[,-1], type="vote")
pc2 <- predict(svm.digits, test[,-1])
table(pc2)
classError(pc2,test[,1]) ## 0.036

###kNN

library(class)

cl <- train[,1]
ptm3 <- proc.time()
knn.digits <- knn(train[,-1], test[,-1], cl, k=5, prob=T, use.all=T)
proc.time() - ptm3 ## (0.19  0.00  1.4)
attributes(knn.digits)

table(knn.digits)
classError(knn.digits, test[,1]) ## 0.0036

## the result using component data is similar to the original,
## but the time used is less than the original one  


#-------------------------------------------------------------------

plotlinearsvm2D=function(svp,xtrain)
## Pretty plot a linear SVM with decision boundary ##
## xtrain should be a 2-dimensional data.
{
	# Define the range of the plot
	# First column is plotted vertically
	yr <- c(min(xtrain[,1]), max(xtrain[,1]))
	# Second column is plotted horizontally
	xr <- c(min(xtrain[,2]), max(xtrain[,2]))

	# Plot the points of xtrain with different signs for positive/negative and SV/non SV
	plot(xr,yr,type='n')
	ymat <- ymatrix(svp)
	points(xtrain[-SVindex(svp),2], xtrain[-SVindex(svp),1], pch = ifelse(ymat[-SVindex(svp)] < 0, 2, 1))
	points(xtrain[SVindex(svp),2], xtrain[SVindex(svp),1], pch = ifelse(ymat[SVindex(svp)] < 0, 17, 16))
	
	# Extract w and b from the model	
	w <- colSums(coef(svp)[[1]] * xtrain[SVindex(svp),])
	b <- b(svp)
	
	# Draw the lines 
	abline(b/w[1],-w[2]/w[1])
	abline((b+1)/w[1],-w[2]/w[1],lty=2)
	abline((b-1)/w[1],-w[2]/w[1],lty=2)
}


cvpred.ksvm <- function(x,y,folds=3,predtype="response",...)
## Return a vector of predictions by cross-validation
## 'predtype' should be one of response (by default), decision or probabilities, depending the prediction we want (SVM label, score or probability, see predict.ksvm())
## Additional parameters are passed to ksvm() to train the SVM
{
 	n <- length(y)
 	ypred <- numeric(n)
 	s <- cv.folds(n,folds)
 	for (i in seq(folds)) {
 		m <- ksvm(x[-s[[i]],],y[-s[[i]]],...)
 		ypred[s[[i]]] <- predict(m,x[s[[i]],],type=predtype)
 		}
 	invisible(ypred)
}



