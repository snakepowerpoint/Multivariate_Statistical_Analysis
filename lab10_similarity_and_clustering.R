setwd('D:/MyGitHub/Multivariate_Statistical_Analysis/data')

#####exercise1-------------------------------------------------
install.packages("HSAUR")
library(HSAUR)
library(xtable)

y <- pottery
yscaled <- t(scale(t(y)))
c <- cor(t(yscaled), method="spearman")
d.1 <- as.dist(1-c)
d.2 <- dist(yscaled, method="euclidean")

hr.11 <- hclust(d.1, method = "single", members=NULL)
hr.12 <- hclust(d.1, method = "complete", members=NULL)
hr.13 <- hclust(d.1, method = "average", members=NULL)

hr.21 <- hclust(d.2, method = "single", members=NULL)
hr.22 <- hclust(d.2, method = "complete", members=NULL)
hr.23 <- hclust(d.2, method = "average", members=NULL)

par(mfrow=c(1,2))
plot(hr.11, hang=-1)
plot(hr.21, hang=-1)

windows()
par(mfrow=c(1,2))
plot(hr.12, hang=-1)
plot(hr.22, hang=-1)

windows()
par(mfrow=c(1,2))
plot(hr.13, hang=-1)
plot(hr.23, hang=-1)

install.packages("fossil")
library(fossil)

cut.11 <- cutree(hr.11, k=3)
cut.21 <- cutree(hr.21, k=3)
cut.12 <- cutree(hr.12, k=4)
cut.22 <- cutree(hr.22, k=4)
cut.13 <- cutree(hr.13, k=4)
cut.23 <- cutree(hr.23, k=3)

## inter-dist
rand.index(cut.11, cut.21)
rand.index(cut.12, cut.22)
rand.index(cut.13, cut.23)
## different distances result in a little difference results

## inter-method based on d.1
rand.index(cut.11, cut.12)
rand.index(cut.11, cut.13)
rand.index(cut.12, cut.13) 
## complete and average are similar under d.1

## inter-method based on d.2
rand.index(cut.21, cut.22)
rand.index(cut.21, cut.23) 
rand.index(cut.22, cut.23)
## three methods based on Euclidean distance give the similar results 

windows()
plot(hr.13, hang=-1)
rect.hclust(hr.13, k=5)

#####exercise2-------------------------------------------------
## we use K-medoids and determine the number k
## distance was calculated with Euclidean
library(cluster)

hp <- pam(d.2, k=3, diss=T)
clusplot(y, hp$clustering, color=TRUE, shade=TRUE, labels=2, lines=0)
plot(hp)

## silihoutte
list.k <- 2:8
avg.sil <- rep(0, time=length(list.k))
for (i in (1:length(list.k))) {
 avg.sil[i] <- pam(d.2, k=list.k[i], diss=T)$silinfo$avg.width
}
print(cbind(list.k, avg.sil)) ## 3 or 4
xtable(cbind(list.k, avg.sil))

## CH index
install.packages("clusterCrit")
library(clusterCrit)

ch <- numeric(0)
for(i in 1:length(list.k)){
 cl <- pam(d.2, k=list.k[i], diss=T)
 ch <- c(ch, intCriteria(as.matrix(yscaled), cl$cluster, 
      "Calinski_Harabasz")$calinski_harabasz)
}
plot(list.k, ch, type="b") ## 4 & 8

## Gap

gap <- clusGap(yscaled, pam, B=100, K.max=8)
plot(gap)
out.gap <- gap$Tab
out.gap
maxSE(out.gap[,3], out.gap[,4], method="globalmax") ## 8
maxSE(out.gap[,3], out.gap[,4], method="Tibs2001SEmax") ## 4

#####exercise3-------------------------------------------------
## By the results above, the number of cluster group is 3,4,8
## under K-medoids cluster method
## now we use K-means cluster method

d.3 <- as.dist(1-cor(t(yscaled))) ## correlation matrix of scaled data


## K-means
## silihoutte
avg.sil.1 <- rep(0, time=length(list.k))
for (i in (1:length(list.k))) {
avg.sil.1[i] <- summary(silhouette(kmeans(yscaled, centers=list.k[i], 
iter.max=200, nstart=25)$cluster, d.3))$avg.width
}
print(cbind(list.k, avg.sil.1)) ## 3
xtable(cbind(list.k, avg.sil.1))

## CH index
ch.1 <- numeric(0)
for(i in 1:length(list.k)){
 cl <- kmeans(yscaled, centers=list.k[i], iter.max=200, nstart=25)
 ch.1 <- c(ch.1, intCriteria(as.matrix(yscaled), cl$cluster, 
      "Calinski_Harabasz")$calinski_harabasz)
}
plot(list.k, ch.1, type="b") ## 4 or 8

## Gap

gap.1 <- clusGap(yscaled, kmeans, B=500, K.max=8)
plot(gap.1)
out.gap.1 <- gap.1$Tab
out.gap.1
maxSE(out.gap.1[,3], out.gap.1[,4], method="globalmax") ## 8
maxSE(out.gap.1[,3], out.gap.1[,4], method="Tibs2001SEmax") ## 3

## the better number of k is 3,4,8 under K-means cluster method


## Spectral clustering
library(kernlab)

hs.3 <- specc(as.matrix(yscaled), centers=3)
hs.4 <- specc(as.matrix(yscaled), centers=4)

hp.3 <- pam(d.3, k=3, diss=T)
hp.4 <- pam(d.3, k=4, diss=T)
hp.8 <- pam(d.3, k=8, diss=T)

hk.3 <- kmeans(d.3, centers=3, iter.max=200, nstart=10)
hk.4 <- kmeans(d.3, centers=4, iter.max=200, nstart=10)
hk.8 <- kmeans(d.3, centers=8, iter.max=200, nstart=10)

rand.index(hs.3, hs.4) ## similar

rand.index(hp.3$clust, hp.4$clust) ## similar
rand.index(hp.3$clust, hp.8$clust)
rand.index(hp.4$clust, hp.8$clust)

rand.index(hk.3$clust, hk.4$clust) ## similar
rand.index(hk.3$clust, hk.8$clust)
rand.index(hk.4$clust, hk.8$clust)

## According to silihoutte value, we can compare K-means with K-medoids
print(cbind(list.k, avg.sil.1)) ## K-means
print(cbind(list.k, avg.sil)) ## K-medoids

## so we use K-means, and choose k=4
clusplot(yscaled, hk.4$cluster, color=TRUE, shade=TRUE, labels=2, lines=0, 
main="K-means")

clusplot(y, hk.4$cluster, color=TRUE, shade=TRUE, labels=2, lines=0, 
main="K-means")


##############################################################
#####exercise1
install.packages("HSAUR")
library(HSAUR)
data(pottery)
?pottery

y <- pottery
## centralize or scale
# scale(t(y)); yscaled <- t(scale(t(y))); apply(yscaled, 1, sd)

# calculate distance
d <- dist(y, method = "euclidean")
# rows of y are observations, columns of y are variables 

# c <- cor(t(y), method="spearman"); d <- as.dist(1-c);
# In order to get the distance based on correlation, we first calculate
# correlation and then turn it into distance. Note that the "cor" function
# calculate the correlations between columns, we need to transpose it 

# hierarchical clustering method
hr <- hclust(d, method = "complete", members=NULL)
plot(hr, hang = -1)



#####exercise2,3
# hierarchical clustering method
hr <- hclust(d, method = "complete", members=NULL)
hc <- cutree(hr, k=3)
plot(hr)
rect.hclust(hr, k=3)
library(cluster)
clusplot(y, hc, color=TRUE, shade=TRUE, labels=2, lines=0, 
main="Complete Linkage")

# K-means
hk <- kmeans(d, centers=3, iter.max=100, nstart=10)
windows()
clusplot(y, hk$cluster, color=TRUE, shade=TRUE, labels=2, lines=0, 
main="K-means")

# K-medoids
hp <- pam(d, k=3, diss=T)
windows()
clusplot(y, hp$clustering, color=TRUE, shade=TRUE, labels=2, lines=0, 
main="K-medoid")

# Spectral clustering
library(kernlab)
hs <- specc(as.matrix(y), centers=3)
windows()
clusplot(y, hs,color=TRUE, shade=TRUE, labels=2, lines=0, main="Spectral")

# compare these results, what do you get?
# Answer: K-means & K-medoids & spectral share the same results



