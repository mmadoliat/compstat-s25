library(bmp); library(jpeg)

pic <- read.bmp("chess.bmp")
svd.pic <- svd(pic)
par(mfrow=c(2,2))
ts.plot(svd.pic$u[,1:2],col=1:2, main="U[,1:2]")
ts.plot(svd.pic$v[,1:2],col=1:2, main="V[,1:2]")
plot(svd.pic$d^2/sum(svd.pic$d^2), main="Scree Plot", type="b")
hist(as.vector(pic))

pic.compress <-  svd.pic$u[,1:2] %*% diag(svd.pic$d[1:2]+10^-8) %*% t(svd.pic$v[,1:2])
writeJPEG(pic.compress,"chess-comp.jpg")

noisy.pic <- pic + rnorm(prod(dim(pic)))
svd.npic <- svd(noisy.pic)
plot(svd.npic$d^2/sum(svd.npic$d^2), main="Scree Plot", log="y", type="b")
pic.compress <-  svd.npic$u[,1:2] %*% diag(svd.npic$d[1:2]) %*% t(svd.npic$v[,1:2])
writeJPEG(noisy.pic,"chess-noisy.jpg")
writeJPEG(pic.compress,"chess-comp.jpg")


pic <- readJPEG('MU.jpeg')
dim(pic)
r <- pic[,,1]
g <- pic[,,2]
b <- pic[,,3]
pic.r.svd <- svd(r)
pic.g.svd <- svd(g)
pic.b.svd <- svd(b)
plot(pic.r.svd$d^2/sum(pic.r.svd$d^2), main="Scree Plot", log="y", type="b")

rgb.svds <- list(pic.r.svd, pic.g.svd, pic.b.svd)

for (j in floor(seq(3, min(dim(pic)[1:2]), length.out = 8))) {
  a <- sapply(rgb.svds, function(i) {
    pic.compress <- i$u[,1:j] %*% diag(i$d[1:j]) %*% t(i$v[,1:j])
  }, simplify = 'array')
  writeJPEG(a, paste('compressed/pic_compressed', '_svd_rank_', round(j,0), '.jpg', sep=''))
}
