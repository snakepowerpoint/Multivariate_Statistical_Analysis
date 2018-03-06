setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####exercise1

t806 <- read.table("T8-6.dat")
names(t806)
names(t806) <- c("Country","100","200","400","800.m","1500.m","5000.m",
"10000.m","Marathon.m")

install.packages("psych")
library(psych)

cov(t806[,2:9])
cortest.bartlett(t806[,2:9])
fa.parallel(t806[,2:9])

t806.fa.1 <- factanal(t806[,2:9], factors=1, scores="regression", 
rotation="varimax", method="mle")
t806.fa.1

library(xtable)
xtable(matrix(t806.fa.1$loadings,1,8))

install.packages("rela")
library(rela)

t806.fa.2 <- paf(as.matrix(t806[,2:9]))
summary(t806.fa.2) ## suggest one factor

xtable(matrix(t806.fa.2$Factor,1,8))

install.packages("GPArotation")
library(GPArotation)

t806.fa.3 <- fa(t806[,2:9], nfactors=1, scores="regression", 
rotate="varimax", fm="mle")
t806.fa.3

xtable(matrix(t806.fa.3$loadings,1,8))

attributes(t806.fa.2)
loadings(t806.fa.1)
t806.fa.2$Factor
loadings(t806.fa.3)
## three types of factor analysis give very similar results,
## every variable has a large loading in factor one
## we can guess that factor one means the ability of athletics

## now we try the case of 3 factors

t806.fa.4 <- factanal(t806[,2:9], factors=3, scores="regression",
rotation="varimax", method="mle")
t806.fa.4
## factor 1 has a large loading in long distance items,
## we say this factor means long distance run endurance
## factor 2 has a large loading in short distance items,
## we say this factor means the ability of short distance run
## factor 3 is difficult to explain

xtable(t(matrix(t806.fa.4$loadings,8,3)))

library(rgl) 

t806.fa.4$scores
plot3d(t806.fa.4$scores, type="n")
text3d(t806.fa.4$scores[,1], t806.fa.4$scores[,2], t806.fa.4$scores[,3], 
texts=1:nrow(t806))

snapshot3d("123.png","png")
## or
rgl.snapshot("123.png", fmt = "png") 

with(t806, plot3d(t806.fa.4$scores, type="s")
identify3d(t806.fa.4$scores[,1], t806.fa.4$scores[,2], 
t806.fa.4$scores[,3],row.names(t806))


#####exercise2

install.packages("plsgenomics")
library(plsgenomics)

data(leukemia)
names(leukemia)
dim(leukemia$gene.names)
dim(leukemia$Y)
dim(leukemia$X) ## what we want to analysis

X <- leukemia$X

X.pc <- prcomp(X, scale.=T) ## use matrix R instead of S
summary(X.pc)

plot(1:38, X.pc$sdev^2, type="b", xlab="Components", ylab="eigenvalue")
plot(X.pc, type="lines") ## choose the first 10 components

dim(X.pc$rotation) ## eigen vectors
head(X.pc$rotation)
dim(X.pc$x) ## x represents values of the first 38 components

attributes(X.pc)
X.pc$sdev

X.L <- X.pc$rotation[,1:10]
X.L <- X.L%*%diag(X.pc$sdev[1:10])
dim(X.L)
head(X.L)

X.L.r <- varimax(X.L)$loadings
attributes(varimax(X.L))

X.U <- diag(diag(cor(X) - X.L.r%*%t(X.L.r)))

inv.cor <- solve(cor(X), tol=1e-25)

X.score <- t(X.L.r)%*%inv.cor%*%t(scale(X)) 
dim(X.score) 
## column means f_j, j in 1 to 38
## row means factor, total factors are 10

plot(t(X.score)[,1:2], type="n", xlab="F1 scores", ylab="F2 scores")
text(t(X.score)[,1], t(X.score)[,2], labels=leukemia$Y, col=leukemia$Y)

par(mfrow=c(3,3))
F <- c()
for(i in 2:10) { 
  nam <- paste("F",i, sep="")
  F[i] <- nam
}

for (i in 2:10){
plot(t(X.score)[,c(1,i)], type="n", xlab="F1", ylab=F[i])
text(t(X.score)[,1], t(X.score)[,i], labels=leukemia$Y, col=leukemia$Y)
i=i+1
}
## two kinds of patients separate

#---------------------useless--------------------------------------
## useless because there exist errors while estimating solve(cor(X))

X.fa.1 <- factanal(X, factors=10, scores="regression", 
rotation="varimax", method="mle")

X.fa.2 <- fa(X, nfactors=10, scores="regression", rotate="varimax", fm="mle")

#---------------------useless--------------------------------------


#####exercise1

library(foreign)
y <- read.spss("http://www.subjectpool.com/ed_teach/y3method/factorexdata05.sav")
x <- as.data.frame(y)
dim(x)
# The data x consists of 538 cases with 102 variables.

for (i in 1:length(x)) { x[,i] <- ifelse(x[,i]==999,NA,x[,i]) }

# The data \verb!x! consists of 538 cases with 102 variables.
# it can be saved as "factorexdata05.txt" by the following line
# write.table(x,"factorexdata05.txt",quote=FALSE,sep="\t",row.names=FALSE)
# if so, the data can be read by:
# x <- read.delim("factorexdata05.txt")

Ps <- x[,4:43] # Extract variables p1-p40
Ps <- subset(Ps, complete.cases(Ps)) 
dim(Ps) # Omit missings (511 cases remain)
res0 <- factanal(Ps,factors=10,scores="reg",rotation="varimax")

install.packages("rela")
library(rela)
res <- paf(as.matrix(Ps))
summary(res)
# Automatically calculate KMO with MSA, determine the number of factors,
# calculate chi-square of Bartlett's sphericity test, communalities and
# factor loadings. Communalities are 1 minus uniquenesses.

attributes(res)
barplot(res$Eigenvalues[,1]) # First column of eigenvalues.
resv <- varimax(res$Factor.Loadings)
# Varimax rotation is possible later.
print(resv)
barplot(sort(colSums(loadings(resv)^2),decreasing=TRUE))
# screeplot using rotated SS loadings.

library(psych)
install.packages("GPArotation")

cortest.bartlett(Ps) # Bartlett's sphericity test.
res2 <- fa.parallel(Ps)
res3 <- fa(Ps, fm="minres", nfactors=8, rotate="varimax")
print(res3)


#####exercise2

install.packages("plsgenomics")
library(plsgenomics)
data(leukemia)
?leukemia
names(leukemia)

#-------------Appendix----------------------------------------------

john = matrix(rnorm(15),ncol=3)
for(i in 1:5) { 
  nam <- paste("b",i, sep="")
  assign(nam, john[i,])
  b[i] <- nam
}


