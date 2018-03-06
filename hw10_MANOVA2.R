setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####7.24
###(a)
###(i)

library(xtable)
library(car)

t110 <- read.table("T1-10.dat", col.names=c("Breed","SalePr","YrHgt",
"FtFrBody","PrctFFB","Frame","BkFat","SaleHt","SaleWt"))
attach(t110)

library(car)

ox.lm <- lm(SaleHt ~ YrHgt + FtFrBody, data=t110)
summary(ox.lm)
xtable(summary(ox.lm))

par(mfrow=c(2,2))
plot(ox.lm)
## mean of residual is close to 0
## normal distribution

str(ox.lm)
plot(ox.lm$residuals)

###(ii)

z0 <- matrix(c(1,50.5,970),3,1)
oxbeta <- ox.lm$coeff
oxre <- ox.lm$residuals
oxz <- cbind(rep(1,76),YrHgt,FtFrBody)

oxlow <- t(z0)%*%oxbeta - qt(1-0.05/2,76-(2+1))*sqrt((1+
t(z0)%*%solve(t(oxz)%*%oxz)%*%z0)*var(oxre)*(76-1)/(76-(2+1)))

oxup <- t(z0)%*%oxbeta + qt(1-0.05/2,76-(2+1))*sqrt((1+
t(z0)%*%solve(t(oxz)%*%oxz)%*%z0)*var(oxre)*(76-1)/(76-(2+1)))

c(oxlow,oxup)

###(b)
###(i)

ox1.lm <- lm(cbind(SaleHt,SaleWt) ~ YrHgt + FtFrBody, data=t110)
summary(ox1.lm)

## three methods to check object
str(ox1.lm)
attributes(ox1.lm)
names(ox1.lm)

ox.manova <- Manova(ox1.lm)
summary(ox.manova)
## SSP means var-cov matrix under hypothesis
## SSPE means var-cov matrix without any condition
## refer to p.307 in textbook
## so we use both YrHgt and FtFrBody

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

xtable.Anova.mlm(ox.manova)

###(ii)

attributes(ox1.lm)
cen <- c(ox1.lm$coeff[,1]%*%z0,ox1.lm$coeff[,2]%*%z0)
shape <- 75*cov(ox1.lm$residuals)/76
r <- c(sqrt((1+t(z0)%*%solve(t(oxz)%*%oxz)%*%z0)*(76/(76-(2+1)))*
(2*(76-(2+1))/(76-(2+2)))*qf(0.95,2,76-(2+2))))

plot(ellipse(cen,shape,r,center.pch=19,center.cex=1.5,draw=FALSE),type='l')
points(x=cen[1],y=cen[2],pch=19,col='red')
lines(c(oxlow,oxup),c(1250,1250),lty=2,col=2,lwd=2)
## prediction ellipse is wider than prediction interval

#####7.26
###(a)
###(i)

t707 <- read.table("T7-7.dat", col.names=c("BL","EM","SF","BS","AFL",
"LFF","FFF","ZST"))

summary(lm(BL ~ AFL + LFF + FFF + ZST, data=t707))
summary(lm(BL ~ LFF + FFF + ZST, data=t707))
BL.lm <- lm(BL ~ LFF + FFF + ZST, data=t707)
## Adjusted R-squared: 0.722

summary(lm(EM ~ AFL + LFF + FFF + ZST, data=t707))
summary(lm(EM ~ AFL + ZST, data=t707))
EM.lm <- lm(EM ~ AFL + ZST, data=t707)
## Adjusted R-squared: 0.7573

summary(lm(SF ~ AFL + LFF + FFF + ZST, data=t707))
summary(lm(SF ~ LFF + FFF + ZST, data=t707))
SF.lm <- lm(SF ~ LFF + FFF + ZST, data=t707)
## Adjusted R-squared: 0.7968

summary(lm(BS ~ AFL + LFF + FFF + ZST, data=t707))
summary(lm(BS ~ LFF + FFF + ZST, data=t707))
BS.lm <- lm(BS ~ LFF + FFF + ZST, data=t707)
## Adjusted R-squared: 0.7444

xtable(BL.lm)
xtable(EM.lm)
xtable(SF.lm)
xtable(BS.lm)

###(ii)

## search outliers by H matrix
attach(t707)
paper.z <- cbind(rep(1,62),AFL,LFF,FFF,ZST)
paper.h <- paper.z%*%solve(t(paper.z)%*%paper.z)%*%t(paper.z)

summary(diag(paper.h))
var(diag(paper.h))
hist(diag(paper.h),freq=F)
which(diag(paper.h)>(0.0806+4*0.0073))
which(diag(paper.h)>0.2) ## observations 46,57,58,60,61

## search outliers by looking error
plot(BL.lm) ## observations 51,52,56
plot(EM.lm)
plot(SF.lm) ## observations 52,56
plot(BS.lm) ## observations 51,52,56

###(iii)

paper.z0 <- matrix(c(1,0.330,45.500,20.375,1.010),5,1)
SF.beta <- lm(SF ~ AFL + LFF + FFF + ZST, data=t707)$coeff
SF.resi <- lm(SF ~ AFL + LFF + FFF + ZST, data=t707)$resi

SF.low <- t(paper.z0)%*%SF.beta - sqrt(61*var(SF.resi)/(62-(4+1)))*sqrt(
1+t(paper.z0)%*%solve(t(paper.z)%*%paper.z)%*%paper.z0)*qt(1-0.05/2,62-(4+1))
SF.up <- t(paper.z0)%*%SF.beta + sqrt(61*var(SF.resi)/(62-(4+1)))*sqrt(
1+t(paper.z0)%*%solve(t(paper.z)%*%paper.z)%*%paper.z0)*qt(1-0.05/2,62-(4+1))
c(SF.low,SF.up)

###(b)
###(i)

paper.lm <- lm(cbind(BL,EM,SF,BS) ~ AFL + LFF + FFF + ZST, data=t707)
paper.lm$coeff
Sigma <- t(paper.lm$resi)%*%paper.lm$resi/62 ## maximum likelihood estimator
Sigma

xtable(paper.lm$coeff)
xtable(Sigma)

###(ii)

## H matrix is the same with (a)(ii),
## so observations 46,57,58,60,61 are outliers
## Also, since multivariate multiple regression model is combination of 
## single-response regression models. So the analysis of residual is the same
## We claim observations 51,52,56 are outliers

###(iii)

Manova(paper.lm)
summary(Manova(paper.lm))
## z1 should be cancled, but we still use all the regressors

Beta <- paper.lm$coeff

pre <- matrix(0,4,2)
for (i in 1:4){
  pre[i,1] <- t(paper.z0)%*%Beta[,i] - sqrt((4*(62-(4+1)))*qf(0.95,4,62-(4+4))*
    (1+t(paper.z0)%*%solve(t(paper.z)%*%paper.z)%*%paper.z0)*62*diag(Sigma)[i]/
    ((62-(4+4))*(62-(4+1))))
  pre[i,2] <- t(paper.z0)%*%Beta[,i] + sqrt((4*(62-(4+1)))*qf(0.95,4,62-(4+4))*
    (1+t(paper.z0)%*%solve(t(paper.z)%*%paper.z)%*%paper.z0)*62*diag(Sigma)[i]/
    ((62-(4+4))*(62-(4+1))))
}
pre
c(SF.low,SF.up)
## the third interval is wider than interval in (a)(iii)


#####################################################
a <- c(12,12,4,5,63,63,23)
b <- c(13,15,7,10,73,83,43)
npts   <- length(a)
shape  <- var(cbind(a, b))
center <- c(mean(a),mean(b))
rconf  <- sqrt(2 * (npts-1) * qf(0.95, 2, npts-2)/(npts*(npts-2)))
rpred  <- sqrt(npts+1)*rconf

conf.elip <- ellipse(center, shape, rconf,draw = FALSE)
pred.elip <- ellipse(center, shape, rpred,draw = FALSE)
plot(pred.elip, type='l')
points(a,b)
lines(conf.elip,col="red")
points(x=center[1],y=center[2],pch=19)

