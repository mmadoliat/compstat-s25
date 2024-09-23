##### Condition Number #####
set.seed(70); 
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
hm <- hilbert(11); x <- 1:nrow(hm); xo <- x+rnorm(length(x),0,.01)
y <- solve(hm)%*%x
ye <- solve(hm)%*%xo

########################## 1D-Basis #########################
rm(list=ls()); require("splines")

set.seed(70); x.range <- c(-180,180); n <- 200; basis <- "bs"; df <- 20; model <- 1; 
t.x <- 0.99*seq(x.range[1],x.range[2],length.out=n)/180*pi;  x <- runif(n,x.range[1],x.range[2])/180*pi;

if (model==1) {
  t.y <- pi^2-t.x^2;  y <- pi^2-x^2+rnorm(n,0,2);  labl <- "TRUE Model: y = pi^2 - x^2"
} else {
  t.y <- sin(t.x);  y <- sin(x)+rnorm(n,0,0.2);  labl <- "TRUE Model: y = sin(x)"
}

if (basis=="poly") {
  B <- outer(x,0:(df-1),"^"); tB <- outer(t.x,0:(df-1),"^");
} else {
  B <- bs(x, df = df); tB <- bs(t.x, df = df); 
}

image(tB); kappa((t(B)%*%B),exact=TRUE)
est <- solve(t(B)%*%B)%*%t(B)%*%y

par(mfrow=c(1,2),mgp=c(1.5,.5,0),mar=c(3,3,1.5,.5),cex=1)
plot(t.x*180/pi,t.y,col=2,ylim=range(y),xlab="x",ylab="y",main=labl,type="l");
points(x*180/pi,y); points((x*180/pi)[order(x)],(B%*%est)[order(x)],type="l",col=3)
ts.plot(tB,col=1:ncol(tB))

