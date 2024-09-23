# Generate a smoothing matrix
get.pen <- function(td) {
  m = length(td);
  h = td[2:m] - td[1:(m-1)]; 
  Q = matrix(0, m, m-1);
  R = matrix(0, m-1, m-1);
  for(k in 2:(m-1))
  {
    Q[k-1,k] = 1/h[k-1];
    Q[k,k] = -1/h[k-1] - 1/h[k];     
    Q[k+1,k] = 1/h[k]
  }
  for(j in 2:(m-2))
  {
    R[j,j] = 1/3 * (h[j-1] + h[j]);
    R[j,j+1] = 1/6 * h[j];
    R[j+1,j] = 1/6 * h[j]
  }
  R[m-1,m-1] = 1/3 * (h[m-2] + h[m-1]);
  s <- solve(R[2:(m-1), 2:(m-1)]) %*% t(Q[1:m, 2:(m-1)]);
  OMEGA = Q[1:m, 2:(m-1)] %*% s;
  EIG.O <- eigen(OMEGA);
  return(list(OMEGA=OMEGA, EIG.O=EIG.O, s=s))
}

# Generate independent of parameter quantities
get.once <- function(EIG.O,X) {if (memory.size()>700) gc();
  once = list(GAMMA=EIG.O$vectors, LAMBDA=diag(EIG.O$values), X.GAMMA=X%*%EIG.O$vectors);
  once
}

# Generate candidates for alphas (page 13 of Dr. Huang's)
get.Alphas <- function(n=101,a=2,s=-75) {return(a^seq(s,s+n))}

# Degrees of freedom in the denominator of GCV
get.df <- function(LAMBDA, alphas=get.Alphas()) {
  i <- 0; df <- array(0, length(alphas));
  for(alpha in alphas)  { i <- i + 1; df[i] =  sum(1/(1+alpha*diag(LAMBDA))) }
  df
}

# Return the GCV for given data matrix, smoothing matrix, list of alphas for k SVD component together
get.GCVk <- function(X, GAMMA, LAMBDA, k=1, alphas=get.Alphas(), trunc.SVD=FALSE) {
  i <- 0; pow <- 0.5; n <- nrow(X); m <- ncol(X); 
  S_p <- array(0, c(m,m,length(alphas))); X.S_p <- array(0, c(n,m,length(alphas)));
  U <- array(0, c(n,k,length(alphas))); V <- array(0, c(m,k,length(alphas)));
  D <- array(0, c(k,k,length(alphas))); GCV <- array(0, length(alphas));
  for(alpha in alphas) {if (memory.size()>700) gc();
	  i <- i + 1; incProgress(1/length(alphas), detail = paste0("2^",log2(alpha)));
   	S_p[,,i] <- GAMMA%*%diag((1/(1+alpha*diag(LAMBDA)))^pow)%*%t(GAMMA);
	  if (trunc.SVD) {SVD.alpha <- propack.svd(X%*%S_p[,,i],k)} else {SVD.alpha <- svd(X%*%S_p[,,i],nu=k,nv=k)}
	  U[,,i] <- SVD.alpha$u[,1:k];
	  if (k==1) D[,,i] <- SVD.alpha$d[1:k] else D[,,i] <- diag(SVD.alpha$d[1:k]);
	  V[,,i] <- S_p[,,i]%*%SVD.alpha$v[,1:k];
	  if (k==1) GCV[i] <- sum((V[,,i]*D[,,i]-t(X)%*%U[,,i])^2) else GCV[i] <- sum((V[,,i]%*%D[,,i]-t(X)%*%U[,,i])^2)
  }
  DF <- get.df(LAMBDA, alphas); GCV <- (1/m)*GCV/((1-(1/m)*DF)^2);
  GCVk = list(U=U, D=D, V=V, GCVs=GCV, DFs=DF, lambdas=alphas);
  GCVk
}

# Find the first k, functional principal component simultaneously
get.FPCk <- function(EIG.O, X, k=1, alphas=get.Alphas(), trunc.SVD=FALSE) {
  once <- get.once(EIG.O, X);
  FPCk <- get.GCVk(X, once$GAMMA, once$LAMBDA, k, alphas, trunc.SVD);
  FPCk
}
