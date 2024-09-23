#######################################################
###       Statistical Computing with R  2e          ###
###       Maria L. Rizzo                            ###
###       Chapman & Hall/CRC The R Series           ###
###       ISBN 9781466553323 - CAT# K15269          ###
###       January 24, 2023                          ###
###                                                 ###
###       R code for Chapter 2                      ###
###       Probability and Statistics Review         ###
#######################################################

### Equations (2.1) and (2.2)

N <- 100; m.X <- v.X <- rep(N)
n <- 1000;
mu <- 5; sig <- 2; df <- 3

for (i in 1:N) {
  Y <- rnorm(n, mean = mu, sd = sig)
  Z <- rchisq(n, df)
  X <- rnorm(n, mean = Y, sd = sqrt(Z))
  m.X[i] <- mean(X);
  v.X[i] <- var(X);
}

mean(m.X);
mean(v.X)

### (2.6) Statistics

N <- 1000; m.X <- v.X <- v.Xb <- rep(N)
n <- 10;
mu <- 5; sig <- 2; df <- 3

for (i in 1:N) {
  Y <- rnorm(n, mean = mu, sd = sig)
  Z <- rchisq(n, df)
  X <- rnorm(n, mean = Y, sd = sqrt(Z))
  m.X[i] <- mean(X);
  v.X[i] <- var(X);
  v.Xb[i] <- sum((X-m.X[i])^2)/n
}

mean(m.X)
mean(v.X)
mean(v.Xb)