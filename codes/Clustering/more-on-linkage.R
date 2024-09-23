set.seed(0)
x = rbind(matrix(rnorm(2*20,mean=-1),ncol=2),
matrix(rnorm(2*30,mean=1),ncol=2))
x = rbind(x,matrix(runif(2*10,min(x),max(x)),ncol=2))
d = dist(x)

# Centroid linkage
tree.cent = hclust(d,method="centroid")
labs.cent = cutree(tree.cent,k=3)

# Minimax linkage
library(protoclust)
tree.mm = protoclust(d)
a = protocut(tree.mm,h=2.5)
labs.mm = a$cl
protos = a$protos

par(mfrow=c(2,2))

plot(tree.cent,labels=FALSE,hang=-1e-10)
abline(h=1.25,lty=2)
plot(tree.mm,labels=rep("",60),hang=-1e-10)
abline(h=2.5,lty=2)

plot(x,col=labs.cent+1)
plot(x,col=labs.mm+1)
points(x[protos,],pch=19,cex=2,col=2:4)

