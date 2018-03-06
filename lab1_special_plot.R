setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

t16<-read.table("T1-6.dat")
library(lattice)

xyplot(V2 ~ V4, data = t16,
groups = V6,
type = c("p", "smooth"), span=.75,
auto.key =list(title = "Iris Data",
x = .15, y=.85, corner = c(0,1),
border = TRUE, lines = TRUE))

xyplot(V3 ~ V5, data = t16,
groups = V6,
type = c("p", "smooth"), span=.75,
auto.key =list(title = "Iris Data",
x = .15, y=.85, corner = c(0,1),
border = TRUE, lines = TRUE))

t19<-read.table("T1-9.dat")
library(rgl)
open3d()
x <- t19$V6
y <- t19$V7
z <- t19$V8
plot3d(x, y, z, col=rainbow(1000))
with(t19, plot3d(V6,V7,V8,type='s', size=2))
identify3d(t19$V6,t19$V7,t19$V8,row.names(t19)) #find 46,11,40

points3d(t19[46,],t19[40,],t19[11,], col="blue")

us<-read.csv("USairpollution.csv")
library(TeachingDemos)
faces2(us[,2:8],labels=as.character(us$X))

library(tourr)
library(andrews)
source("starcoord.R")
source("mmnorm.R")
source("circledraw.R")
source("radviz2d.R")
andrews(flea, type = 4, clr = 5, ymax = 2, main = "Type = 4")
parallelplot(~flea[1:6], flea, groups = species,
horizontal.axis = FALSE, scales = list(x = list(rot = 90)))
starcoord(flea,class = T)
radviz2d(flea)







