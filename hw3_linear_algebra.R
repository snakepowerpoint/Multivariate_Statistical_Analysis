A<-array(c(1,2,2,1,-2,2),dim=c(3,2))
AT<-t(A)
B<-AT%*%A
eigen(B)

E<-A%*%AT
eigen(E)

svd(A)
svd(A)$u %*% diag(svd(A)$d) %*% t(svd(A)$v) # check the answer


library(x)
x=xtable(A) 
print(x, floating=FALSE, tabular.environment="bmatrix",
hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)