setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####exercise1
###(1)

tbody <- read.table("body.dat", col.names=c("ID","weight","height",
"neck","chest","hip","thigh","biceps","wrist","agegp"))
attach(tbody)

n1 <- table(agegp)[1]
n2 <- table(agegp)[2]

factor <- factor(tbody$agegp,levels=c("young","old"),labels=c('1','2'))
body <- tbody[-c(1,10)]
body <- cbind(body,factor)

x1 <- body[factor==1, -9]
x2 <- body[factor==2, -9]
p <- ncol(x1)

m1 <- apply(x1, 2, mean)
m2 <- apply(x2, 2, mean)
Sp <- ((n1-1)*var(x1)+(n2-1)*var(x2))/(n1+n2-2)

T2 <- ((n1*n2)/(n1+n2))* (t(m1-m2) %*% solve(Sp) %*% (m1-m2) )
Fstat <- ((n1+n2-p-1)*T2)/((n1+n2-2)*p)
pvalue <- 1-pf(Fstat,p, n1+n2-p-1)

print(paste("Hotelling T^2 =", round(T2,4),
"F=", round(Fstat,4), "P-value =", round(pvalue,4) ))
## reject H_0
## check every variable

l <- rep(0,8)
u <- rep(0,8)
for (i in 1:8){
l[i] = m1[i]-m2[i]-sqrt((n1+n2)*diag(Sp)[i]/(n1*n2))*qt(1-0.05/(2*p),n1+n2-2)
u[i] = m1[i]-m2[i]+sqrt((n1+n2)*diag(Sp)[i]/(n1*n2))*qt(1-0.05/(2*p),n1+n2-2)
}
cbind(l,u)
## so variable 2 and 4 are different in mean

library(xtable)
xtable(cbind(l,u))

###(2)
###normality

library(MASS)
library(energy)

## H_0 is normal assumption
mvnorm.etest(x1) ## can not reject normal assumption
mvnorm.etest(x2) ## can not reject normal assumption

###homogeneity of covariance

## H_0 is homogeneity of covariance

Boxm <- function(X, cl, alpha=0.05) { 

if (alpha <= 0 || alpha >= 1) 
  stop('significance level must be between 0 and 1') 
g = nlevels(cl) ## Number of groups. 
n = table(cl) ## Vector of groups-size. 
N = nrow(X) 
p = ncol(X) 
bandera = 2 
if (any(n >= 20))
  bandera = 1 
## Partition of the group covariance matrices. 

covList <- tapply(as.matrix(X), rep(cl, ncol(X)), function(x, nc) cov(matrix(x, nc = nc)),
                  ncol(X))
deno = sum(n) - g 
suma = array(0, dim=dim(covList[[1]])) 
for (k in 1:g) 
  suma = suma + (n[k] - 1) * covList[[k]] 
Sp = suma / deno ## Pooled covariance matrix. 
Falta=0 
for (k in 1:g) 
  Falta = Falta + ((n[k] - 1) * log(det(covList[[k]]))) 

MB = (sum(n) - g) * log(det(Sp)) - Falta ## Box's M statistic. 
suma1 = sum(1 / (n[1:g] - 1)) 
suma2 = sum(1 / ((n[1:g] - 1)^2)) 
C = (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) * (g - 1))) * 
  (suma1 - (1 / deno)) ## Computing of correction factor. 
if (bandera == 1)
  { 
    X2 = MB * (1 - C) ## Chi-square approximation. 
    v = as.integer((p * (p + 1) * (g - 1)) / 2) ## Degrees of freedom. 
    ## Significance value associated to the observed Chi-square statistic. 
    P = pchisq(X2, v, lower=FALSE)  #RM: corrected to be the upper tail 
 
    if (P >= alpha) { 
 
    } else { 
 
    } 
    return(list(MBox=MB, ChiSq=X2, df=v, pValue=P)) 
  }
else
  { 
    ## To obtain the F approximation we first define Co, which combined to 
    ## the before C value are used to estimate the denominator degrees of 
    ## freedom (v2); resulting two possible cases. 
    Co = (((p-1) * (p+2)) / (6 * (g-1))) * (suma2 - (1 / (deno^2))) 
    if (Co - (C^2) >= 0) { 
      v1 = as.integer((p * (p + 1) * (g - 1)) / 2) ## Numerator DF. 
      v21 = as.integer(trunc((v1 + 2) / (Co - (C^2)))) ## Denominator DF. 
      F1 = MB * ((1 - C - (v1 / v21)) / v1) ## F approximation. 
      ## Significance value associated to the observed F statistic. 
      P1 = pf(F1, v1, v21, lower=FALSE) 

      if (P1 >= alpha) { 
         
      } else { 
         
      } 
      return(list(MBox=MB, F=F1, df1=v1, df2=v21, pValue=P1)) 
    } else { 
      v1 = as.integer((p * (p + 1) * (g - 1)) / 2) ## Numerator df. 
      v22 = as.integer(trunc((v1 + 2) / ((C^2) - Co))) ## Denominator df. 
      b = v22 / (1 - C - (2 / v22)) 
      F2 = (v22 * MB) / (v1 * (b - MB)) ## F approximation. 
      ## Significance value associated to the observed F statistic. 
      P2 = pf(F2, v1, v22, lower=FALSE) 
         
      if (P2 >= alpha) { 
        cat('Covariance matrices are not significantly different.\n') 
      } else { 
        cat('Covariance matrices are significantly different.\n') 
      } 
      return(list(MBox=MB, F=F2, df1=v1, df2=v22, pValue=P2)) 
    } 
  }
}

X <- rbind(x1,x2)
group <- factor(rep(1:2,c(n1,n2)))
Boxm(X, group, 0.05) 
## reject homogeneity of covariance
## so their covariance are different
## however, we look their cov

cov(x1)
cov(x2) ## their covariances are actually familiar

xtable(cov(x1))
xtable(cov(x2))

-------------------------------------------------------
#####exercise2
###(1)

rmvnormmix<-function (n, p, mus, D){
mus <- matrix(mus, nrow = nrow(mus))
r <- NCOL(mus)
z <- sample(1:2, n, replace = TRUE, prob =c(p,1-p))
if(p==0) rand<-mvrnorm(n,mu=mus[2,],D)
else if (p==1) rand<-mvrnorm(n,mu=mus[1,],diag(r))
else{
rand<-mvrnorm(n, mu =rep(0,r), diag(r))+mus[z, ]
rand[z==2,]<-rand[z==2,]%*%sqrt(D)
}
rand
}

N1 <- 100
N2 <- 100 #sample size
q <- 10 #dimension of variables
Mu1 <- rep(0,q)
Mu2 <- rep(5,q)

set.seed(1234)
Y1<-rmvnormmix(N1,p=1,mus=rbind(Mu1,Mu2),2*diag(q))
Y2<-rmvnormmix(N2,p=0,mus=rbind(Mu1,Mu2),2*diag(q))

##normality test
mvnorm.etest(Y1) ## reject normality
mvnorm.etest(Y2) ## reject mormality

#Box M test of homogeneity of covariance
Y<-rbind(Y1,Y2)
group<-factor(rep(1:2,c(N1,N2)))
Boxm(Y,group,0.05)

#actual type I error, nominal level is 0.05
idx<-NA
for( i in 1:1000){
Y1<-rmvnormmix(N1,p=0.5,mus=rbind(Mu1,Mu2),2*diag(q))
Y2<-rmvnormmix(N2,p=0.5,mus=rbind(Mu1,Mu2),2*diag(q))
Y<-rbind(Y1,Y2)
group<-factor(rep(1:2,c(N1,N2)))
idx<-c(idx,Boxm(Y,group,0.05)$pValue)
}

idx<-idx[-1]
sum(idx<0.05)/1000 ## is not familiar to 0.05

###(2)

boxnorm <- function(n1,n2,p1,p2,q){
 mu1<-rep(0,q)
 mu2<-rep(5,q)
 idx<-NA
  for(i in 1:1000){
   Y1<-rmvnormmix(n1,p1,mus=rbind(mu1,mu2),2*diag(q))
   Y2<-rmvnormmix(n2,p2,mus=rbind(mu1,mu2),2*diag(q))
   Y<-rbind(Y1,Y2)
   group<-factor(rep(1:2,c(n1,n2)))
   idx<-c(idx,Boxm(Y,group,0.05)$pValue) 
  }
 idx<-idx[-1]
 sum(idx<0.05)/1000
}

boxnorm(100,100,0.2,0.8,10)
boxnorm(100,100,0,1,10)
boxnorm(100,100,1,0,10)
boxnorm(100,100,1,1,10)
boxnorm(100,100,0,0,10)

##################################################
##################################################
##################################################
###exercise1

steel <- read.table("steel.dat", header=F, 
col.names=c("temperature", "yield", "strength"))
attach(steel)

# now we perform the two-sample Hotelling T^2-test
# compute step-by-step
#----------------------------------------------------
# number of observations for each group:
n1<-table(temperature)[1]
n2<-table(temperature)[2]

# Splitting the data matrix (while removing the 1th column)
# into two subsets, one for group 1 and another for group 2:
x1 <- steel[temperature==1, -1]
x2 <- steel[temperature==2, -1]
p <- ncol(x1)

# Sample mean vectors for each group:
m1 <- apply(x1, 2, mean)
m2 <- apply(x2, 2, mean)

# "pooled" sample covariance matrix:
Sp <- ((n1-1)*var(x1)+(n2-1)*var(x2))/(n1+n2-2)

# Hotelling T^2, the F-statistic, and the P-value:
T2 <- ((n1*n2)/(n1+n2))* (t(m1-m2) %*% solve(Sp) %*% (m1-m2) )
Fstat <- ((n1+n2-p-1)*T2)/((n1+n2-2)*p)
pvalue <- 1-pf(Fstat,p, n1+n2-p-1)

print(paste("Hotelling T^2 =", round(T2,4),
"F=", round(Fstat,4), "P-value =", round(pvalue,4) ))

#use HotellingsT2 in ICSNP package
#----------------------------------------------------
library(ICSNP)
HotellingsT2(steel[temperature == 1, -1], steel[temperature == 2, -1])

#alternatively,
#----------------------------------------------------
summary(manova(cbind(yield, strength) ~ temperature),test="Hotelling")

