setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####exercise1
###(1)

###normality

iris

library(energy)
mvnorm.etest(iris[1:50,1:4]) ## non-normal
mvnorm.etest(iris[51:100,1:4]) ## cannot reject normal
mvnorm.etest(iris[101:150,1:4]) ## cannot reject normal

###var-cov matrix

cov(iris[1:50,1:4])
cov(iris[51:100,1:4])
cov(iris[101:150,1:4]) ## seems like dissimilar

BoxMTest <- function(X, cl, alpha=0.05) { 
if (alpha <= 0 || alpha >= 1) 
  stop('significance level must be between 0 and 1') 
g = nlevels(cl) ## Number of groups. 
n = table(cl) ## Vector of groups-size. 
N = nrow(X) 
p = ncol(X) 
bandera = 2 
if (any(n >= 20))
  bandera = 1 
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
    cat('------------------------------------------------\n'); 
    cat(' MBox Chi-sqr. df P\n') 
    cat('------------------------------------------------\n') 
    cat(sprintf("%10.4f%11.4f%12.i%13.4f\n", MB, X2, v, P)) 
    cat('------------------------------------------------\n') 
    if (P >= alpha) { 
      cat('Covariance matrices are not significantly different.\n') 
    } else { 
      cat('Covariance matrices are significantly different.\n') 
    } 
    return(list(MBox=MB, ChiSq=X2, df=v, pValue=P)) 
  }
else
  { 
    Co = (((p-1) * (p+2)) / (6 * (g-1))) * (suma2 - (1 / (deno^2))) 
    if (Co - (C^2) >= 0) { 
      v1 = as.integer((p * (p + 1) * (g - 1)) / 2) ## Numerator DF. 
      v21 = as.integer(trunc((v1 + 2) / (Co - (C^2)))) ## Denominator DF. 
      F1 = MB * ((1 - C - (v1 / v21)) / v1) ## F approximation. 
      ## Significance value associated to the observed F statistic. 
      P1 = pf(F1, v1, v21, lower=FALSE) 
      cat('\n------------------------------------------------------------\n') 
      cat(' MBox F df1 df2 P\n') 
      cat('------------------------------------------------------------\n') 
      cat(sprintf("%10.4f%11.4f%11.i%14.i%13.4f\n", MB, F1, v1, v21, P1)) 
      cat('------------------------------------------------------------\n') 
      if (P1 >= alpha) { 
        cat('Covariance matrices are not significantly different.\n') 
      } else { 
        cat('Covariance matrices are significantly different.\n') 
      } 
      return(list(MBox=MB, F=F1, df1=v1, df2=v21, pValue=P1)) 
    } else { 
      v1 = as.integer((p * (p + 1) * (g - 1)) / 2) ## Numerator df. 
      v22 = as.integer(trunc((v1 + 2) / ((C^2) - Co))) ## Denominator df. 
      b = v22 / (1 - C - (2 / v22)) 
      F2 = (v22 * MB) / (v1 * (b - MB)) ## F approximation. 
      ## Significance value associated to the observed F statistic. 
      P2 = pf(F2, v1, v22, lower=FALSE) 
      
      cat('\n------------------------------------------------------------\n') 
      cat(' MBox F df1 df2 P\n') 
      cat('------------------------------------------------------------\n') 
      cat(sprintf('%10.4f%11.4f%11.i%14.i%13.4f\n', MB, F2, v1, v22, P2)) 
      cat('------------------------------------------------------------\n') 
      
      if (P2 >= alpha) { 
        cat('Covariance matrices are not significantly different.\n') 
      } else { 
        cat('Covariance matrices are significantly different.\n') 
      } 
      return(list(MBox=MB, F=F2, df1=v1, df2=v22, pValue=P2)) 
    } 
  }
}

factor <- factor(iris[,5],levels=c("setosa","versicolor","virginica"),
labels=c('1','2','3'))

BoxMTest(iris[,1:4],factor,0.05) ## var-cov matrices are different
BoxMTest(iris[,1:4],iris[,5],0.05) 
## another way(without generate factor)

###(2)

trteff.diff<-function (data, level=0.95) {
# Bonferroni Based Simultaneous CI for Treatments?? Difference for one-way MANOVA
#data: Nx(p+1) numeric matrix or dataframe of data, N is the total sample size, p  is number of variables,
# the last column is the factors
#level: confidence level of interval,default=0.95

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
cat("\n Bonferroni Based Simultaneous CI for Treatments?? Difference \n")
SCI = data.frame(Trt = C, Estimate = round(tau,3), LowerCI = round(lower, 3), UpperCI = round(upper,3))
labs<-NA
for(i in 1:p)
 for(k in 1:(g-1))
   for(l in (k+1):g)
     labs<-c(labs,paste("tau",k,"-",l,i,sep=""))
rownames( SCI)<-labs[-1]
print(SCI)
}

trteff.diff(iris)

library(xtable)
xtable(trteff.diff(iris))

#####exercise2
###(a)

cm560 <- c(10.35,13.41,7.78,10.40,17.78,10.40)
cm720 <- c(25.93,38.63,25.15,24.25,41.45,29.20)
breed <- c(5,6,8,5,6,8)
nutrient <- c(1,1,1,2,2,2)

t632 <- data.frame(cbind(cm560,cm720,breed,nutrient))
t632$breed <- as.factor(t632$breed)
t632$nutrient <- as.factor(t632$nutrient)

#----------------------------------------------------
tree.lmm <- lm(cbind(cm560,cm720) ~ breed*nutrient, data=t632)
tree.lmm
summary(tree.lmm)
attributes(tree.lmm)
#----------------------------------------------------

tree.lm <- lm(cbind(cm560,cm720) ~ breed + nutrient, data=t632)
summary(tree.lm)

library(car)
tree.manova <- Manova(tree.lm)
tree.manova
summary(tree.manova)

###(b)

tree.lm1 <- lm(cm560 ~ breed + nutrient, data=t632)
tree.anova1 <- Anova(tree.lm1)
tree.anova1

tree.lm2 <- lm(cm720 ~ breed + nutrient, data=t632)
tree.anova2 <- Anova(tree.lm2)
tree.anova2

#####exercise3

t633 <- read.table("T6-18.dat", col.names=c("cm560", "cm720", "breed",
"time", "copy"))
t633 <- t633[-5]
t633$breed <- as.factor(t633$breed)
t633$time <- as.factor(t633$time)

lm633.1 <- lm(cbind(cm560,cm720) ~ breed*time, data=t633)
summary(lm633.1)

lm633.2 <- lm(cbind(cm560,cm720) ~ time*breed, data=t633)
## different turn

## typeI SS
manova(lm633.1)
manova(lm633.2) ## they are different

summary(manova(lm633.1))
summary(manova(lm633.2)) ## the same result

## typeII SS
Manova(lm633.1)
summary(Manova(lm633.1))

## typeIII SS
Manova(lm(cbind(cm560,cm720) ~ breed*time, data=t633,
contrasts=list(breed=contr.sum, time=contr.sum)), type="III")

Manova(lm633.1, type="III") ## without contrast is different
summary(Manova(lm633.1, type="III"))

#------------------------------------------------------------

tree.lm1 <- lm(cm560 ~ breed + nutrient, data=t632)
john <- Anova(tree.lm1)
john
tree.anova1

xtable.Anova.mlm <- function (x, ...) {
  test <- x$test
  repeated <- x$repeated
  ntests <- length(x$terms)
  tests <- matrix(NA, ntests, 4)
  if (!repeated)
    SSPE.qr <- qr(x$SSPE)
  for (term in 1:ntests) {
    eigs <- Re(eigen(qr.coef(if (repeated) qr(x$SSPE[[term]]) else
      SSPE.qr,
      x$SSP[[term]]), symmetric = FALSE)$values)
    tests[term, 1:4] <- switch(test, Pillai = stats:::Pillai(eigs,
      x$df[term], x$error.df), Wilks = stats:::Wilks(eigs,
        x$df[term], x$error.df), `Hotelling-Lawley` = stats:::HL(eigs,
          x$df[term], x$error.df), Roy = stats:::Roy(eigs,
            x$df[term], x$error.df))
  }
  ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
  ok <- !is.na(ok) & ok
  tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3],
    tests[ok, 4], lower.tail = FALSE))
  rownames(tests) <- x$terms
  colnames(tests) <- c("Df", "test stat", "approx F", "num Df",
    "den Df", "Pr(>F)")
  tests <- structure(as.data.frame(tests), heading = paste("\nType ",
    x$type, if (repeated)
      " Repeated Measures", " MANOVA Tests: ", test, " test
statistic",
    sep = ""), class = c("anova", "data.frame"))
  #    print(tests)
  #    invisible(x)
  xtable(tests)
}


