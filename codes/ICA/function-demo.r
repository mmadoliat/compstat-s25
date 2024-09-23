library("fastICA"); ncomp=4

sensors <- as.matrix(read.csv("data.csv", header=FALSE));
data <- t(sensors); data <- scale(data);

ts.plot(data)
ts.plot(t(data))

ica <- fastICA(data,ncomp)
pca <- princomp(data,ncomp, corr=FALSE)
SVD <- svd(data)

plot(pca$loadings[,1],pca$loadings[,2])

par(mfrow=c(2,2))
for(i in 1:ncomp) plot(ica$S[,i],type="l")
for(i in 1:ncomp) plot(ica$A[i,],type="l")
for(i in 1:ncomp) plot(SVD$u[,i],type="l")
for(i in 1:ncomp) plot(SVD$v[,i],type="l")
