##########################Regular grid, tensor spline#########################
rm(list=ls()); require("splines"); library("scatterplot3d")

set.seed(70); x.range <- y.range <- c(-180,180); n <- 20; basis <- "bs"; df <- 8; grid <- "regular"; model <- 1;
t.x <- seq(x.range[1],x.range[2],length.out=n)/180*pi;  t.y <- seq(y.range[1],y.range[2],length.out=n)/180*pi;

if (grid=="regular") {
  t.data <- cbind(rep(t.x,each=length(t.y)),rep(t.y,length(t.x))); x <- t.data[,1]; y <- t.data[,2];
} else {
  # Does not work with tensor product
  x <- runif(n^2,x.range[1],x.range[2])/180*pi; y <- runif(n^2,y.range[1],y.range[2])/180*pi; 
}

if (model==1) {
	t.z <- matrix(pi^2-t.data[,1]^2-t.data[,2]^2,nr=n);  z <- pi^2-x^2-y^2+rnorm(n^2,0,2);  labl <- "TRUE Model: z = pi^2 - x^2 - y^2"
} else {
	t.z <- matrix(sin(t.data[,2]),nr=n);  z <- sin(y)+rnorm(n^2,0,0.2);  labl <- "TRUE Model: z = sin(y)"
}
data <- cbind(x,y); 
scatterplot3d(x,y,z)

if (basis=="poly") {
  B <- outer(t.x,0:(df-1),"^"); tB <- outer(t.x,0:(df-1),"^");
} else {
  B <- bs(t.x, df = df); tB <- bs(t.x, df = df); 
}
B <- B%x%B; tB <- tB%x%tB

image(tB); kappa((t(B)%*%B),exact=TRUE)
for (i in 2:ncol(B)) persp(matrix(B[,i],nr=length(t.x)))
est <- solve(t(B)%*%B)%*%t(B)%*%z

#pdf(file=paste(model+2,"test.pdf",sep=""),width=12,height=6);
par(mfrow=c(1,2),mgp=c(1.5,.5,0),mar=c(3,3,1.5,.5),cex=1)
 persp(t.x*180/pi,t.y*180/pi,t.z,theta=-40,phi=40,zlim=range(z),xlab="x",ylab="y",zlab="z",main=labl);
 contour(t.x*180/pi,t.y*180/pi,t.z,xlab="x",ylab="y",col=2)

 persp(t.x*180/pi,t.y*180/pi,matrix(B%*%est,nr=length(t.x)),theta=-40,phi=40,zlim=range(z),xlab="x",ylab="y",zlab="z",main="\n LSE: Based on 400 samples + Gaussian noise");
 contour(t.x*180/pi,t.y*180/pi,matrix(B%*%est,nr=length(t.x)),xlab="x",ylab="y",col=2)
#dev.off()

################################# Triangulation ##################################
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
rm(list=ls()); gc(); library(INLA); source("BBP.R")

### Obtaining Triangles ###
x.range <- y.range <- c(-180,180)
loc.bnd = matrix(c(x.range[1],y.range[1], x.range[2],y.range[1], x.range[2],y.range[2], x.range[1],y.range[2]), 4, 2, byrow=TRUE);
segm.bnd = inla.mesh.segment(loc.bnd)
ell = max(x.range[2]-x.range[1],y.range[2]-y.range[1])
mes <- round(c(sqrt(2*ell^2)+.01,ell+.01,sqrt(2*(ell/2)^2)+.01,(ell/2)+.01,sqrt(2*(ell/4)^2)+.01,(ell/4)+.01,sqrt(2*(ell/8)^2)+.01,(ell/8)+.01,sqrt(2*(ell/16)^2)+.01),3); mesh <- list();
for (i in 1:length(mes)) {mesh[[i]] = inla.mesh.create(boundary=segm.bnd,refine=list(max.edge=mes[i]))}
par(mfrow=c(3,3),mgp=c(-1,0,0),mar=c(1,0,3,0),cex=.75)
    for (i in 1:9) {plot(mesh[[i]]); mtext(paste(mesh[[i]]$n,"vertices,",nrow(mesh[[i]]$graph$tv),"triangles for max.edge",mes[i]))}

set.seed(70); d <- 3; r <- 1; n <- 400; mesh1 <- mesh[[2]]; mesh2 <- mesh[[3]]; model <- 1
t.x <- 0.99*seq(x.range[1],x.range[2],16)/180*pi; t.y <- seq(y.range[1],y.range[2],16)/180*pi;  
t.data <- cbind(rep(t.x,each=length(t.y)),rep(t.y,length(t.x)))
x <- runif(n,x.range[1],x.range[2])/180*pi; y <- runif(n,y.range[1],y.range[2])/180*pi;

if (model==1) {
	t.z <- matrix(pi^2-t.data[,1]^2-t.data[,2]^2,nr=length(t.x));  z <- pi^2-x^2-y^2+rnorm(n,0,2);  labl <- "TRUE Model: z = pi^2 - x^2 - y^2"
} else {
	t.z <- matrix(sin(t.data[,2]),nr=length(t.x));  z <- sin(y)+rnorm(n,0,0.2);  labl <- "TRUE Model: z = sin(y)"
}
 data <- cbind(x,y); plot(data)
 scatterplot3d(x,y,z)

tB1 <- basis(mesh1,d,t.data)
tB2 <- basis(mesh2,d,t.data)

B <- basis(mesh1,d,data)
est <- (solve(t(B)%*%B)%*%t(B)%*%z)

B <- basis(mesh2,d,data)
  H <- smoothness(mesh2,d,r,type="regular",intercept=FALSE)
	qrdecom=qr(t(H));	Q=qr.Q(qrdecom,TRUE); C=Q[,(qrdecom$rank+1):ncol(Q)]
	BC <- B%*%C; tBc <- tB2%*%C;
 estC <- (solve(t(BC)%*%BC)%*%t(BC)%*%z)

  aH <- smoothness(mesh2,d,r,type="angular",intercept=FALSE)
	qrdecom=qr(t(aH)); Q=qr.Q(qrdecom,TRUE); A=Q[,(qrdecom$rank+1):ncol(Q)]
	BA <- B%*%A; tBa <- tB2%*%A;
 estA <- (solve(t(BA)%*%BA)%*%t(BA)%*%z)

z.lim <- range(z); nrz <- nrow(z); ncz <- ncol(z); 

#pdf(file=paste(model,"test.pdf",sep=""),width=12,height=6);
par(mfrow=c(1,2),mgp=c(1.5,.5,0),mar=c(3,3,1.5,.5),cex=1)
 persp(t.x*180/pi,t.y*180/pi,t.z,theta=-40,phi=40,zlim=z.lim,xlab="x",ylab="y",zlab="z",main=labl);
 contour(t.x*180/pi,t.y*180/pi,t.z,xlab="x",ylab="y",col=2,asp=1,xaxt="n",yaxt="n")
 plot(mesh1,add=T); axis(c(1,2), c(-180,-90,0,90,180), round(c(-pi,-pi/2,0,pi/2,pi),2)); axis(2, c(-180,-90,0,90,180), round(c(-pi,-pi/2,0,pi/2,pi),2))

 persp(t.x*180/pi,t.y*180/pi,matrix(tB1%*%est,nr=length(t.x)),theta=-40,phi=40,zlim=z.lim,xlab="x",ylab="y",zlab="z",main="\n LSE: No Boundary Conditions\n Based on 100 samples + Gaussian noise");
 contour(t.x*180/pi,t.y*180/pi,matrix(tB1%*%est,nr=length(t.x)),xlab="x",ylab="y",col=2,asp=1,xaxt="n",yaxt="n")
 plot(mesh1,add=T); axis(c(1,2), c(-180,-90,0,90,180), round(c(-pi,-pi/2,0,pi/2,pi),2)); axis(2, c(-180,-90,0,90,180), round(c(-pi,-pi/2,0,pi/2,pi),2))

 persp(t.x*180/pi,t.y*180/pi,matrix(tBc%*%estC,nr=length(t.x)),theta=-40,phi=40,zlim=z.lim,xlab="x",ylab="y",zlab="z",main="LSE: Interior Boundary Conditions");
 contour(t.x*180/pi,t.y*180/pi,matrix(tBc%*%estC,nr=length(t.x)),xlab="x",ylab="y",col=2,asp=1,xaxt="n",yaxt="n")
 plot(mesh1,add=T); axis(c(1,2), c(-180,-90,0,90,180), round(c(-pi,-pi/2,0,pi/2,pi),2)); axis(2, c(-180,-90,0,90,180), round(c(-pi,-pi/2,0,pi/2,pi),2))

 persp(t.x*180/pi,t.y*180/pi,matrix(tBa%*%estA,nr=length(t.x)),theta=-40,phi=40,zlim=z.lim,xlab="x",ylab="y",zlab="z",main="LSE: Angular-Interior Boundary Conditions");
 contour(t.x*180/pi,t.y*180/pi,matrix(tBa%*%estA,nr=length(t.x)),xlab="x",ylab="y",col=2,asp=1,xaxt="n",yaxt="n")
 plot(mesh1,add=T); axis(c(1,2), c(-180,-90,0,90,180), round(c(-pi,-pi/2,0,pi/2,pi),2)); axis(2, c(-180,-90,0,90,180), round(c(-pi,-pi/2,0,pi/2,pi),2))
#dev.off()

 