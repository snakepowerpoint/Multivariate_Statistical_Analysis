library(MASS)
sigma1 <- matrix(c(2,2^(-1/2),2^(-1/2),1),2,2)
a<-mvrnorm(n=1000,mu=c(0,2),sigma1)
colnames(a)<-c("X1","X2")

install.packages("mixtools")
library(mixtools)
plot(a)
ellipse(mu<-colMeans(a),sigma<-cov(a),alpha=.5,col='red')
