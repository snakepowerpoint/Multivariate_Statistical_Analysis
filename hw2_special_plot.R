setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

library(rgl)
x<-read.table("T11-4.dat")
colnames(x)<-c("CF.TD","NI.TA","CA.CL","CA.NS","Type")
colors <- ifelse(x$Type==0, 'red','blue') 
with(x, plot3d(CF.TD, NI.TA, CA.CL, col=colors, type='s', size=2))

y=x[c(1,2,3,4)]
cor(y)

m<-read.table("T6-10.dat")
colnames(m)<-c("X1","X2","X3","Type")
colors1 <- ifelse(m$Type=='gasoline', 'red','blue')
with(m, plot3d(X1, X2, X3, col=colors1, type='s', size=2))

n<-read.table("T11-9.dat")
colnames(n)<-c("brand","firm","cal","pro","fat","na","fiber",
"carb","s","k","group")
stars(n)

install.packages("aplpack")
library(aplpack)

install.packages("symbols")
library(symbols)
symbol(n,type="face")

z<-read.table("T12-4.dat")
colnames(z)<-c("X1","X2","X3","X4","X5","X6","X7","X8","company")
stars(z)

a<-read.table("T1-10.dat")
colnames(a)<-c("Breed","salePr","YrHgt","YrFrBody","PrctFFB","Frame",
"BkFat","SaleHt","SaleWt")
apply(a,2,mean)
cov(a)
cor(a)

with(a, plot3d(Breed, Frame, BkFat, col=colors1, type='s', size=2))









