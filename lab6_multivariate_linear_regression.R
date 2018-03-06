setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####exercise1
###(1)

library(xtable)

t706 <- read.table("T7-6.dat", col.names=c("TOT","AMI","GEN",
"AMT","PR","DIAP","QRS"))

blue.lm <- lm(cbind(TOT,AMI) ~ GEN + AMT + PR + DIAP + QRS, data=t706)
summary(blue.lm)
xtable(blue.lm$coeff)

###(2)

## multivariate multiple regression
resid(blue.lm)
S <- t(resid(blue.lm))%*%resid(blue.lm)/blue.lm$df

## multivariate regression
blue.lm1 <- lm(TOT ~ GEN + AMT + PR + DIAP + QRS, data=t706)
t(resid(blue.lm1))%*%resid(blue.lm1)/blue.lm1$df

blue.lm2 <- lm(AMI ~ GEN + AMT + PR + DIAP + QRS, data=t706)
t(resid(blue.lm2))%*%resid(blue.lm2)/blue.lm2$df

t(resid(blue.lm1))%*%resid(blue.lm2)/blue.lm2$df
S ## exactly the same

xtable(S)

#####exercise2
###(1)

library(car)

manova.blue <- Manova(blue.lm)
manova.blue ## PR, DIAP, QRS is not significant
summary(manova.blue)
xtable(manova.blue)

###(2)

manova.blue ## gender is a significant variable

blue.lm3 <- lm(cbind(TOT,AMI) ~ GEN, data=t706)
linearHypothesis(blue.lm3, matrix(c(0,1),1,2))
linearHypothesis(blue.lm3, matrix(c(1,0),1,2))
linearHypothesis(blue.lm, "GEN0 = GEN1")

#-----------------------------------------------------------
#####exercise1

attach(iris)
names(iris)
names(iris) <- c("SL", "SW", "PL", "PW", "SPP")

mod.iris <- lm(cbind(SL, SW, PL, PW) ~ SPP, data=iris)
mod.iris

summary(mod.iris)
Sigma<-t(resid(mod.iris))%*%resid(mod.iris)/mod.iris$df

#for purpose of comparison
mod1<-lm(SL~SPP,data=iris)
summary(mod1)
t(resid(mod1))%*%resid(mod1)/mod1$df

mod2<-lm(SW~SPP,data=iris)
summary(mod2)
t(resid(mod2))%*%resid(mod2)/mod2$df

t(resid(mod1))%*%resid(mod2)/mod2$df

#####exercise2

library(car)
manova.iris <- Anova(mod.iris)
summary(manova.iris)

#test for differences between
#setosa and the average of versicolor and virginica

linearHypothesis(mod.iris, "0.5*SPPversicolor + 0.5*SPPvirginica")
#or
linearHypothesis(mod.iris, matrix(c(0,0.5,0.5),nrow=1,ncol=3))

#test for differences between versicolor and virginica:
linearHypothesis(mod.iris, "SPPversicolor = SPPvirginica")
#or
linearHypothesis(mod.iris, matrix(c(0,1,-1),nrow=1,ncol=3))

#----------------------------------------------------------------

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