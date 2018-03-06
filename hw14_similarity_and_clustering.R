setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####12.5
###(a)

m1205 <- matrix(c(0,1,11,5,1,0,2,3,11,2,0,4,5,3,4,0),4,4)
as.dist(m1205)

plot(hclust(as.dist(m1205), method='single'))
plot(hclust(as.dist(m1205), method='complete'))
plot(hclust(as.dist(m1205), method='average'))

library(xtable)
xtable(m1205,digits=0)

#####12.16
###(a)

t109 <- read.table("T1-9.dat")
names(t109)
names(t109) <- c("Country","100","200","400","800.m","1500.m","3000.m", 
"Marathon.m")
row.names(t109) <- t109[,1]
t109 <- t109[,-1]

dist.109 <- dist(t109, method="euclidian")
xtable(as.matrix(dist.109))

###(b)

## single linkage
single.109 <- hclust(dist.109, method='single')
plot(single.109, labels=row.names(t109), ylab="Distance")

windows() # opening new window while keeping previous one open

## complete linkage
complete.109 <- hclust(dist.109, method='complete')
plot(complete.109, labels=row.names(t109), ylab="Distance")

## more clearer graph
single.109 <- as.dendrogram (single.109)
plot(single.109)

complete.109 <- as.dendrogram (complete.109)
plot(complete.109)

## or
install.packages('ape')
library(ape)

single.109 <- hclust(dist.109, method='single')
complete.109 <- hclust(dist.109, method='complete')

plot(as.phylo(single.109), cex = 0.9)
plot(as.phylo(single.109), type = "cladogram", cex = 0.9)

###(c)

install.packages("clusterCrit")
library(clusterCrit)

## determine k by CH index
list.k <- 2:10
ch <- numeric(0)
for(k in list.k){
 cl <- kmeans(t109, centers=k, iter.max=100, nstart=25)
 ch <- c(ch,intCriteria(as.matrix(t109),cl$cluster,
        "Calinski_Harabasz")$calinski_harabasz)
}
plot(list.k, ch, type="l")
which.max(ch)+1 ## seems to be 6

library(cluster)

## For what value of k chosen by silhouette method?
a <- c()
for (i in 2:10){
 john <- summary(silhouette(cutree(single.109,k=i),dist.109))$avg.width
 a <- c(a,john)
}
which.max(a)+1 ## 2

a <- c()
for (i in 2:10){
 john <- summary(silhouette(cutree(complete.109,k=i),dist.109))$avg.width
 a <- c(a,john)
}
which.max(a)+1 ## 2

a <- c()
for (i in 2:10){
 john <- summary(silhouette(kmeans(t109, centers=i, iter.max=100, 
              nstart=25)$cluster, dist.109))$avg.width
 a <- c(a, john)
}
which.max(a)+1 ## 2

## For what value of k in maxSE?
gap.cars <- clusGap(t109, kmeans , B=500, K.max=10)
plot(gap.cars)
out.gap <- gap.cars$Tab
out.gap
maxSE(out.gap[,3],out.gap[,4],method="globalmax")
maxSE(out.gap[,3],out.gap[,4],method="Tibs2001SEmax")

## For what value of k does the elbow of the plot occur?
n <- length(t109[,1])
wss1 <- (n-1)*sum(apply(t109,2,var))
wss <- numeric (0)
for(i in list.k) {
  W <- sum(kmeans(t109,i)$withinss)
  wss <- c(wss,W)
}
wss <- c(wss1,wss)
plot(c(1,list.k),wss,type='l',xlab='Number of clusters', 
ylab='Within-groups sum-of-squares', lwd=2) ## seems to be 4




kmean.109.2 <- kmeans(t109, center=2, iter.max=100, nstart=25)
kmean.109.2

kmean.109.4 <- kmeans(t109, center=4, iter.max=100, nstart=25)
kmean.109.6 <- kmeans(t109, center=6, iter.max=100, nstart=25)

install.packages("fossil")
library(fossil)

rand.index(kmean.109.2$clus,kmean.109.4$clus)
rand.index(kmean.109.2$clus,kmean.109.6$clus)
rand.index(kmean.109.4$clus,kmean.109.6$clus) ## 4 or 6 is ok

## the result is consistent with the linkages procedures

##---------------------------------------------------------------------
##----------------------        Appendix         ----------------------
##---------------------------------------------------------------------


## put histograms on the diagonal
panel.hist  <- function(x, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
     }


## put bivariate density estimation on the upper panels,
panel.kernel.density  <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, ...) 
{
    points(x, y,pch = pch, col =col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
    	xl<-list(x1=x,x2=y)
    	xx<-cbind(x,y)
        d <- bkde2D(xx,bandwidth=sapply(xl,dpik))
     contour(x=d$x1,y=d$x2,d$fhat,add=T,...)  
}
}

library(KernSmooth)


