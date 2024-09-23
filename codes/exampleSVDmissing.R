set.seed(30)
x <- seq(0,2*pi,length.out=100)
cs <- sort(rnorm(20))
X <- true <- cs%*%t(cos(x))
ts.plot(t(true),col=1:length(cs), main="TRUE")

Xo <- X <- X + rnorm(prod(dim(X)),sd=0.05)
ts.plot(t(Xo),col=1:length(cs), main="OBSERVED")

miss.rate <- 0.8
ind <- sample(1:prod(dim(X)), size = prod(c(miss.rate,dim(X))), replace = FALSE)
X[ind] <- NA; XO <- X;
ts.plot(t(XO),col=1:length(cs), main=paste0("Missing : ", miss.rate*100, "%"))

X[ind] <- mean(X,na.rm = T)
ts.plot(t(X),col=1:length(cs))
for (i in 1:2^12) {
  svd.X <- svd(X)
  X.tilde <- svd.X$u[,1]%*%diag(svd.X$d[1],1)%*%t(svd.X$v[,1])
  X[ind] <- X.tilde[ind]
  if (!(log2(i)%%1)) ts.plot(t(X),col=1:length(cs), main=paste("iter",i))
}

filled.contour(t(true),col=1:length(cs), main="TRUE")
filled.contour(t(XO),col=1:length(cs), main=paste0("Missing : ", miss.rate*100, "%"))
filled.contour(t(X),col=1:length(cs), main=paste("iter",i))
filled.contour(t(true),col=1:length(cs), main="TRUE")
filled.contour(t(Xo),col=1:length(cs), main="OBSERVED")

