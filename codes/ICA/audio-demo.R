library("fastICA"); library("tuneR");

orig <- list();
orig[[1]] <- readWave("orig/cult-band.wav")
orig[[2]] <- readWave("orig/siren.wav")
orig[[3]] <- readWave("orig/brandenburg6.wav")
orig[[4]] <- readWave("orig/xmas-song.wav")
orig[[5]] <- readWave("orig/cnn-male.wav")
orig[[6]] <- readWave("orig/cnn-female.wav")
orig[[7]] <- readWave("orig/ica-intro.wav")
orig[[8]] <- readWave("orig/ica-more.wav")
orig[[9]] <- readWave("orig/read-aloud.wav")
smp <- sample(2:length(orig))

o1 <- orig[[smp[1]]]; writeWave(o1, filename = "out/orig1.wav")
o2 <- orig[[smp[2]]]; writeWave(o2, filename = "out/orig2.wav")
o3 <- orig[[smp[3]]]; writeWave(o3, filename = "out/orig3.wav")

play(o1); o1 <- o1@left
play(o2); o2 <- o2@left
play(o3); o3 <- o3@left

a <- runif(3); a <- a/sum(a)
x1 <- a[1]*(o1) + a[2]*(o2) + a[3]*(o3); sens1 <- Wave(x1, samp.rate=8000, bit=8)
x2 <- a[2]*(o1) + a[3]*(o2) + a[1]*(o3); sens2 <- Wave(x2, samp.rate=8000, bit=8)
x3 <- a[3]*(o1) + a[1]*(o2) + a[2]*(o3); sens3 <- Wave(x3, samp.rate=8000, bit=8)
X <- cbind(x1,x2,x3);

play(sens1); writeWave(sens1, filename = "out/sens1.wav")
play(sens2); writeWave(sens2, filename = "out/sens2.wav")
play(sens3); writeWave(sens3, filename = "out/sens3.wav")

pcaS <- princomp(X)$scores; 
pcaS <- pcaS - min(pcaS); pcaS <- min(X) + pcaS*((max(X)-min(X))/max(pcaS))

u1 <- Wave(pcaS[,1], samp.rate=8000, bit=8)
u2 <- Wave(pcaS[,2], samp.rate=8000, bit=8)
u3 <- Wave(pcaS[,3], samp.rate=8000, bit=8)

play(u1); writeWave(u1, filename = "out/pca1.wav")
play(u2); writeWave(u2, filename = "out/pca2.wav")
play(u3); writeWave(u3, filename = "out/pca3.wav")
screeplot(princomp(X))

icaS <- fastICA(X,n.comp=3)$S
icaS <- icaS - min(icaS); icaS <- min(X) + icaS*((max(X)-min(X))/max(icaS))

u1 <- Wave(icaS[,1], samp.rate=8000, bit=8)
u2 <- Wave(icaS[,2], samp.rate=8000, bit=8)
u3 <- Wave(icaS[,3], samp.rate=8000, bit=8)

play(u1); writeWave(u1, filename = "out/ica1.wav")
play(u2); writeWave(u2, filename = "out/ica2.wav")
play(u3); writeWave(u3, filename = "out/ica3.wav")
