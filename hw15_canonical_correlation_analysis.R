setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####10.10
###(a)

R.1010 <- scan()
1.0 0.615 -0.111 -0.266
0.615 1.0 -0.195 -0.085
-0.111 -0.195 1.0 -0.269
-0.266 -0.085 -0.269 1.0

R.1010 <- matrix(R.1010,4,4)
H.1 <- eigen(matrix.msqrt(R.1010[1:2,1:2])%*%R.1010[1:2,3:4]%*%
solve(R.1010[3:4,3:4])%*%R.1010[3:4,1:2]%*%matrix.msqrt(R.1010[1:2,1:2]))

sqrt(H.1$values)

###(b)

## Note that, R.1010[3:4,3:4]^(-1) != solve(R.1010[3:4,3:4])
## refer to text pp.51

H.1$vectors[,1]%*%matrix.msqrt(R.1010[1:2,1:2]) ## U_1

H.2 <- eigen(matrix.msqrt(R.1010[3:4,3:4])%*%R.1010[3:4,1:2]%*%
solve(R.1010[1:2,1:2])%*%R.1010[1:2,3:4]%*%matrix.msqrt(R.1010[3:4,3:4]))

H.2$vectors[,1]%*%matrix.msqrt(R.1010[3:4,3:4]) ## V_1

## X_1 is primary variable in U_1
## V_1 means punishment index
## we see that punishment is correlated with nonprimary homicides
## but not primary homicides


#####10.13
###(a)

R.1 <- scan()
1.0 0 0 0 0
0.754 1.0 0 0 0 
-0.690 -0.712 1.0 0 0
-0.446 -0.515 0.323 1.0 0
0.692 0.412 -0.444 -0.334 1.0

R.1 <- matrix(R.1,5,5)
R.1 <- R.1 + t(R.1) - diag(1,5)

R.2 <- scan()
-0.605 -0.722 0.737 0.527 -0.383
-0.479 -0.419 0.361 0.461 -0.505
0.780 0.542 -0.546 -0.393 0.737
-0.152 -0.102 0.172 -0.019 -0.148

R.2 <- matrix(R.2,5,4)

R.3 <- scan()
1.0 0 0 0
0.251 1.0 0 0
-0.490 -0.434 1.0 0
0.250 -0.079 -0.163 1.0

R.3 <- matrix(R.3,4,4)
R.3 <- R.3 + t(R.3) - diag(1,4)

R <- matrix(0,9,9)
R[1:5,1:5]=R.1
R[1:5,6:9]=R.2
R[6:9,1:5]=t(R.2)
R[6:9,6:9]=R.3

H.1 <- eigen(matrix.msqrt(R.1)%*%R.2%*%solve(R.3)%*%t(R.2)%*%
matrix.msqrt(R.1))
rho <- H.1$values

n <- 138
p <- 5
q <- 4
chi <- numeric()
onepercent <- numeric()
for (i in 1:5){
 chi[i]= -(n-1-(p+q+1)/2)*log(prod((1-rho)[i:5]))
 onepercent[i] <- qchisq(0.99,(p-i+1)*(q-i+1))
}
print(cbind(chi, onepercent))

## the first two canonical correlations are significant

H.2 <- eigen(matrix.msqrt(R.3)%*%t(R.2)%*%solve(R.1)%*%R.2%*%
matrix.msqrt(R.3))

t(H.1$vectors[,1:2])%*%matrix.msqrt(R.1)
t(H.2$vectors[,1:2])%*%matrix.msqrt(R.3)

library(xtable)

xtable(t(H.1$vectors[,1:2])%*%matrix.msqrt(R.1))
xtable(t(H.2$vectors[,1:2])%*%matrix.msqrt(R.3))

###(b)

t(H.1$vectors[,1])%*%matrix.msqrt(R.1)
t(H.2$vectors[,1])%*%matrix.msqrt(R.3)

## U_1 shows the difference between good qualities and bad ones
## V_1 is difficult to explain

###(c)

A <- rbind(t(H.1$vectors[,1])%*%matrix.msqrt(R.1), 
t(H.1$vectors[,2])%*%matrix.msqrt(R.1), 
t(H.1$vectors[,3])%*%matrix.msqrt(R.1), 
t(H.1$vectors[,4])%*%matrix.msqrt(R.1), 
t(H.1$vectors[,5])%*%matrix.msqrt(R.1))

sum(diag(matrix(solve(A)[,1],5,1)%*%solve(A)[,1]))/5
## the contribution of U_1 to X^(1)
## or sum(((A%*%R.1)[1,])^2)/5

B <- rbind(t(H.2$vectors[,1])%*%matrix.msqrt(R.3), 
t(H.2$vectors[,2])%*%matrix.msqrt(R.3), 
t(H.2$vectors[,3])%*%matrix.msqrt(R.3), 
t(H.2$vectors[,4])%*%matrix.msqrt(R.3))

sum(diag(matrix(solve(B)[,1],4,1)%*%solve(B)[,1]))/4
## the contribution of V_1 to X^(2)
## or sum(((B%*%R.3)[1,])^2)/4


#####12.19

t1213 <- scan()
  0.000 0 0 0 0 0 0 0 0
  2.202  0.000 0 0 0 0 0 0 0
  1.004  2.025  0.000 0 0 0 0 0 0 
  1.108  1.943  0.233  0.000 0 0 0 0 0 
  1.122  1.870  0.719  0.541  0.000 0 0 0 0 
  0.914  2.070  0.719  0.679  0.539  0.000 0 0 0 
  0.914  2.186  0.452  0.681  1.102  0.916  0.000 0 0
  2.056  2.055  1.986  1.990  1.963  2.056  2.027  0.000 0
  1.608  1.722  1.358  1.168  0.681  1.005  1.719  1.991  0.000

t1213 <- matrix(t1213,9,9)
t1213 <- t1213 + t(t1213)

fit.1 <- cmdscale(t1213, eig=T, k=2)
attributes(fit.1)

library(MASS)

sts <- rep(0,6)
for (i in 1:6){
 sts[i] <- isoMDS(t1213, k=i)$stress
}
plot(1:6, sts, type="b")

fit.2 <- isoMDS(t1213, k=5)
fit.2

plot(fit.2$points[,1:2], type="n", xlim=1.2*range(fit.2$points[,1]))
text(fit.2$points[,1], fit.2$points[,2], labels=c("Age1","Age8","Age3",
"Age4","Age6","Age5","Age2","Age9","Age7"), col="red")

## This figure shows that the data can be illustrated in low-dimensional
## because the primary difference between the data is time factor.
## The stress figure also coincides the inference above.

john <- princomp(fit.2$points)
plot(john$scores[,1:2], type="n", xlim=1.2*range(john$scores[,1]))
text(john$scores[,1], john$scores[,2], labels=c("Age1","Age8","Age3",
"Age4","Age6","Age5","Age2","Age9","Age7"), col="red")
## exactly the same!!!


#--------------------------------------------------------------
#              Appendix
#--------------------------------------------------------------

matrix.sqrt <- function(x){
 vec <- eigen(x)$vectors
 val <- eigen(x)$values
 h <- matrix(0,nrow(x),nrow(x))
 for (i in 1:nrow(x)){
  ai <- sqrt(val[i])*vec[,i]%*%t(vec[,i])
  h <- h+ai
 }
 return(h)
}

matrix.msqrt <- function(x){
 vec <- eigen(x)$vectors
 val <- eigen(x)$values
 h <- matrix(0,nrow(x),nrow(x))
 for (i in 1:nrow(x)){
  ai <- (sqrt(val[i])^(-1))*vec[,i]%*%t(vec[,i])
  h <- h+ai
 }
 return(h)
}





