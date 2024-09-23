data <- c(rnorm(100,1),rnorm(100,4)); n <- length(data); k <- ans <- 1; B <- 2; eps <- 10^(-5);
par <- list(); par$pi <- par$mu <- par$sigma2 <- list();
par$pi[[1]] <- par$mu[[1]] <- par$sigma2[[k]] <- rep(1/B,B); 
q <- quantile(data,seq(0,1,1/B));
for (j in 1:B) {
  dat <- data[data>=q[j] & data<q[(j+1)]]; 
  par$mu[[k]][j] <- mean(dat); par$sigma2[[k]][j] <- var(dat)
}

p.hat <- matrix(0, nr=n, nc=B);
while (ans==1) {
  for (j in 1:B) {p.hat[,j] <- dnorm(data,mean = par$mu[[k]][j], sd=sqrt(par$sigma2[[k]][j]))}
  p.hat <- p.hat/apply(p.hat,1,sum); k <- k + 1;
  par$pi[[k]] <- apply(p.hat,2,mean);
  par$mu[[k]] <- apply(p.hat*data,2,sum)/(n*par$pi[[k]])
  par$sigma2[[k]] <- apply(p.hat*outer(data,par$mu[[k]],"-")^2,2,sum)/(n*par$pi[[k]])
  if (sum(abs(par$pi[[k]]-par$pi[[(k-1)]])) < eps) ans <- 0
}

Xs <- seq(q[1],q[B+1],length.out=100); d.hat <- matrix(0, nr=length(Xs), nc=B)
for (j in 1:B) {d.hat[,j] <- par$pi[[k]][j]*dnorm(Xs,mean = par$mu[[k]][j], sd=sqrt(par$sigma2[[k]][j]))}
d.hat <- apply(d.hat,1,sum); hist(data, probability = T); points(Xs,d.hat,type="l")
