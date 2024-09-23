#######################################################################     METHOD     ###########################################################################

####################          INITIALIZATION USING KDE         ################
initial=function(Xs,mesh,d=3,r=1,type="angular",intercept=TRUE,dim=length(Xs)){
	X.G <- cbind(rep(Xs[[1]]$x,length(Xs[[1]]$y)),rep(Xs[[1]]$y,each=length(Xs[[1]]$x))); len=length(Xs);
	W <- matrix(0,nr=length(Xs[[1]]$z),nc=len)
	for (i in 1:len) W[,i] <- as.vector(log(Xs[[i]]$z))
	W <- scale(W, center=apply(W,2,mean), scale=FALSE)
	B <- basis(mesh,d,X.G)
	H <- smoothness(mesh,d,r,type=type,intercept)
		qrdecom=qr(t(H))
		Q=qr.Q(qrdecom,TRUE)
		C=Q[,(qrdecom$rank+1):ncol(Q)]
	BC=B%*%C
	smooth <- solve(t(BC)%*%BC)%*%t(BC)
	Theta.tA <- smooth%*%W
	svd.Theta.tA <- svd(Theta.tA);
	k <- min(sum(svd.Theta.tA$d/svd.Theta.tA$d[1]>.Machine$double.eps),dim);
	Theta <- svd.Theta.tA$u[,1:k]; A <- (svd.Theta.tA$v%*%diag(svd.Theta.tA$d))[,1:k];
	A <- A%*%diag(sign(Theta[1,])); Theta <- Theta%*%diag(sign(Theta[1,])); 
	K1=energy(mesh,d,r); K1C=t(C)%*%K1%*%C;
	return(list(W=W, BC=BC, Theta=Theta, A=A, K1C=K1C, H=H, B=Matrix(B), C=C, smooth=smooth, Theta.tA=Theta.tA, hatC=BC%*%smooth, K1=K1, svd.Theta.tA=svd.Theta.tA))
}

#######################         Evaluating the Penalized log-Likelihood     ###############
plogLike=function(Theta=init$Theta,A=init$A,lambda=0) {
	w <- (init$BC%*%Theta%*%t(A))
	if (length(Xs[[1]][[1]][[1]])>1) Xs <- Xs[[1]]
	temp <- (c(Xs[[1]]$x[-1],Xs[[1]]$x[1]+(x.range[2]-x.range[1]))-Xs[[1]]$x)%*%t(c(Xs[[1]]$y[-1],Xs[[1]]$y[1]+(y.range[2]-y.range[1]))-Xs[[1]]$y)
	W <- list(); logLike <- 0; len=length(Xs);
	for (i in 1:len) {
		W[[i]]=once[[i]]$BC%*%Theta%*%A[i,]
		logLike <- logLike + sum(W[[i]])-nrow(W[[i]])*log(sum(exp(w[,i])*temp))
	}
	plogLike <- -2*logLike+lambda*sum(diag(t(Theta)%*%init$K1C%*%Theta))
	return(plogLike)
}

#######################           Iterations for a given tuning parameter          #############
iterate=function(lambda=0) {
	eps <- 0.0001; n.iter <- 1000;
	Theta <- list(); A <- list(); diff <- eps + 1; l <- 1;
	Theta[[l]] <- init$Theta; A[[l]] <- init$A;
	while (abs(diff) > eps && l < n.iter) {
		ll <- 0; repit <- TRUE;
		A.chng <- matrix(0,nrow=nrow(A[[l]]),ncol=ncol(A[[l]])); A[[(l+1)]] <- A.chng;
		Theta.chng <- matrix(0,nrow=nrow(Theta[[l]]),ncol=ncol(Theta[[l]])); Theta[[(l+1)]] <- Theta.chng;
		temp1 <- array(0,dim=c(nrow(init$K1C),nrow(init$K1C),ncol(Theta[[l]]))); tmp1i <- temp1;
		temp2 <- array(0,dim=c(nrow(init$K1C),ncol(Theta[[l]])));
		for (i in 1:nrow(A[[l]])) {
			temp3 <- nrow(once[[i]]$BC)*t(Theta[[l]])%*%once[[i]]$V.B%*%Theta[[l]]
			  if (qr(temp3)$rank<nrow(temp3)) {
				lll <- 0; temp4 <- temp3;
				while (qr(temp4)$rank<nrow(temp4)) {
					temp4 <- temp3 + diag(.Machine$double.eps*(10^lll), nrow(temp3)); lll <- lll + 1
				}
				temp3 <- temp4
			  }
			A.chng[i,] <- solve(temp3)%*%(t(Theta[[l]])%*%once[[i]]$b_E.b);
			for (k in 1:ncol(Theta[[l]])) {
				temp1[,,k] <- temp1[,,k] + nrow(once[[i]]$BC)*A[[l]][i,k]^2*once[[i]]$V.B
				temp2[,k] <- temp2[,k] + A[[l]][i,k]*once[[i]]$b_E.b
				if (i==nrow(A[[l]])) {
					if (qr(temp1[,,k] + lambda*init$K1C)$rank < nrow(temp1[,,k])) {print("choose smaller lambda"); break()}
					tmp1i[,,k] <- solve(temp1[,,k] + lambda*init$K1C)
					Theta.chng[,k] <- tmp1i[,,k]%*%(temp2[,k]-lambda*init$K1C%*%Theta[[l]][,k])
				}
			}
		}
		while (repit) {
			tau <- (0.5)^ll; 
			for (i in 1:nrow(A[[l]])) {A[[(l+1)]][i,] <- A[[l]][i,] + tau*A.chng[i,]}
			for (k in 1:ncol(Theta[[l]])) {Theta[[(l+1)]][,k] <- Theta[[l]][,k] + tau*Theta.chng[,k]}
			if (plogLike(Theta[[(l+1)]],A[[(l+1)]],lambda) <= plogLike(Theta[[l]],A[[l]],lambda)) {repit <- FALSE} else {ll <- ll + 1}
		}
# browser()
		diff <- (plogLike(Theta[[(l+1)]],A[[(l+1)]],lambda) - plogLike(Theta[[l]],A[[l]],lambda));
		print(paste("l=",l," half=",ll," pLogL=",round(plogLike(Theta[[l]],A[[l]],lambda))," diff=",round(diff,-log10(eps)),sep=""));
		if(diff>eps) {browser()};
		l <- l + 1; 
	}
	svd.Theta.tA <- svd(Theta[[(l)]]%*%t(A[[(l)]]))
	Theta[[(l)]] <- svd.Theta.tA$u[,1:ncol(Theta[[l]])]; A[[(l)]] <- (svd.Theta.tA$v%*%diag(svd.Theta.tA$d))[,1:ncol(Theta[[l]])];
	A[[(l)]] <- A[[(l)]]%*%diag(sign(Theta[[(l)]][1,])); Theta[[(l)]] <- Theta[[(l)]]%*%diag(sign(Theta[[(l)]][1,])); 
	return (list(Theta=Theta[[l]],A=A[[l]],lambda=lambda,plogLike=plogLike(Theta[[(l)]],A[[(l)]],lambda),temp1=temp1,tmp1i=tmp1i))
}

############## Akaike information criterion ###########
AIC <- function(lambda=10^-10) {
	ans <- iterate(lambda); df <- 0;
	for (k in 1:ncol(ans$Theta)) {df <- df + sum(diag(ans$tmp1i[,,k]%*%ans$temp1[,,k]))}
	AIC <- plogLike(Theta=ans$Theta,A=ans$A,lambda=0) + 2*df
	return(AIC)
}
#######################################################################     END OF METHOD     ###################################################################



#######################################################################     TRIANGULATION     ###################################################################

################## xy2bary: Transform Cartesian coordinates to Barycentric coordinates. ######
xy2bary=function(mesh.vert,xy){
  if (is.vector(xy) || nrow(xy)==1) xy <- matrix(c(1,xy[1:2]),nc=3,nr=1) else xy <- cbind(1,xy[,1:2]);
  if (is.matrix(mesh.vert)) {bary <- t(solve(t(cbind(1,mesh.vert[,1:2]))) %*% t(xy))} else {
	bary <- matrix(0,nr=nrow(xy),nc=3)
	for (i in 1:nrow(mesh.vert$graph$tv)) {
		vert.i <- mesh.vert$loc[mesh.vert$graph$tv[i,],]
		result.i <- t(solve(t(cbind(1,vert.i[,1:2]))) %*% t(xy))
		indx <- apply(round(result.i,15),1,min)>=0
		bary[indx,] <- result.i[indx,]
	}
  }
  return(bary)
}

################## Bernstein Basis Polynomials ########
basis=function(mesh,d,data){
  n.basis=(d+1)*(d+2)/2; indx.all <- array(0,nrow(data))
	b.matrix=matrix(0,nrow(data),n.basis*nrow(mesh$graph$tv))
  for (i in 1:nrow(mesh$graph$tv)){
    result <- xy2bary(mesh$loc[mesh$graph$tv[i,],],data)
    indx <- (sign(apply(result,1,min))==1)&(!indx.all); indx.all <- indx|indx.all;
    if(sum(indx)){
      for(j in 1:n.basis){
	  c_se=rep(0,n.basis)
	  c_se[j]=1
	  b.matrix[indx,(i-1)*n.basis+j]=loceval(result[indx,1],result[indx,2],result[indx,3],c_se)
	}
    }
  }
  return(b.matrix)
}

################## Used in Basis (local Triangles) #######
loceval=function(lam1_le,lam2_le,lam3_le,bcoef_le){
  nc_le=length(bcoef_le)
  d_le=degree(nc_le)
  for (j in 1:d_le) bcoef_le=de_cast_step(lam1_le,lam2_le,lam3_le,bcoef_le)
  y=t(bcoef_le)
  return(y)
}

################## De Casteljau's Algorithm #######
de_cast_step=function(lam1_bo,lam2_bo,lam3_bo,Bin_bo){
  Bin_bo=as.matrix(Bin_bo)
  m_bo=nrow(Bin_bo)
  d_bo=degree(m_bo)
  n_bo=length(lam1_bo)
  result=indices(d_bo)
  I_bo=result[,1];J_bo=result[,2];K_bo=result[,3]
  result=indices(d_bo-1)
  I1_bo=result[,1];J1_bo=result[,2];K1_bo=result[,3]
  index1_bo=locate(cbind(I1_bo+1,J1_bo,K1_bo),cbind(I_bo,J_bo,K_bo))
  index2_bo=locate(cbind(I1_bo,J1_bo+1,K1_bo),cbind(I_bo,J_bo,K_bo))
  index3_bo=locate(cbind(I1_bo,J1_bo,K1_bo+1),cbind(I_bo,J_bo,K_bo))
  if (ncol(Bin_bo)==1) Bout_bo=Bin_bo[index1_bo]%*%t(lam1_bo)+Bin_bo[index2_bo]%*%t(lam2_bo)+Bin_bo[index3_bo]%*%t(lam3_bo) else {
    if(length(lam1_bo)>1){
	Bout_bo=Bin_bo[index1_bo,]%*%diag(lam1_bo,length(lam1_bo),length(lam1_bo))+Bin_bo[index2_bo,]%*%diag(lam2_bo,length(lam2_bo),
	  length(lam2_bo))+Bin_bo[index3_bo,]%*%diag(lam3_bo,length(lam3_bo),length(lam3_bo))
    } else {
	Bout_bo=Bin_bo[index1_bo,]*lam1_bo+Bin_bo[index2_bo,]*lam2_bo+Bin_bo[index3_bo,]*lam3_bo
    }
  }
  return(Bout_bo) 
}

################# Used in De Casteljau's Algorithm ######
degree=function(m_dg){
  d_dg=(-3+sqrt(8*m_dg+1))/2
  return(d_dg)
}

################# Used in De Casteljau's Algorithm ######
indices=function(d_in){
  m_in=(d_in+1)*(d_in+2)/2
  I_in=rep(0,m_in)
  J_in=I_in
  K_in=I_in
  Mark_in=1
  for (j in seq(d_in,0,-1)){
	I_in[Mark_in:(Mark_in+j)]=seq(j,0,-1)
	J_in[Mark_in:(Mark_in+j)]=0:j
	K_in[Mark_in:(Mark_in+j)]=(d_in-j)*rep(1,j+1)
	Mark_in=Mark_in+j+1
  }
  return(cbind(I_in,J_in,K_in))
}

################# Used in De Casteljau's Algorithm ######
locate=function(matrix1_lc,matrix2_lc){	
  colnames(matrix1_lc)=NULL
  colnames(matrix2_lc)=NULL
  n1_lc=nrow(matrix1_lc)
  n2_lc=nrow(matrix2_lc)
  ind_lc=rep(0,n1_lc)
  for(j in 1:n1_lc){
    for(k in 1:n2_lc){
	if (all(matrix1_lc[j,]==matrix2_lc[k,])) {ind_lc[j]=k;break}
    }
  }
  return(ind_lc)
}
#######################################################################     END OF TRIANGULATION     ###################################################################



#######################################################################         SMOOTHNESS           ###################################################################

################## smoothness: Define smoothness matrix H such that Hc=0. #########
smoothness=function(mesh_sm,d_sm,r_sm,type="angular",intercept=FALSE){
	V_sm <- mesh_sm$loc; T_sm <- mesh_sm$graph$tv;

	# Add artificial triangles in the top and right of the mesh for angular data
	if (type=="angular") {
	  low.tri <- which(T_sm%in%which(V_sm[,2]==-180),arr.ind=TRUE) %% nrow(T_sm); low.tri[low.tri==0] <- nrow(T_sm);
	  left.tri <- which(T_sm%in%which(V_sm[,1]==-180),arr.ind=TRUE) %% nrow(T_sm); left.tri[left.tri==0] <- nrow(T_sm);
	  flag.low <- setequal(V_sm[V_sm[,2]==180,1],V_sm[V_sm[,2]==-180,1]) && sum(V_sm[,2]==180)
	  flag.left <- setequal(V_sm[V_sm[,1]==180,2],V_sm[V_sm[,1]==-180,2]) && sum(V_sm[,1]==180)
	  if (flag.low) {
		low.tri <- low.tri[duplicated(low.tri)]; 
		for (j in 1: length(low.tri)) {
		  V_tmp <- V_sm[T_sm[low.tri[j],],]; V_tmp[,2] <- V_tmp[,2] + 360
		  T_tmp2 <- NULL;
		  for (k in 1:3) {
			v_tmp <- which(apply(t(V_sm)==V_tmp[k,],2,sum)==3)
			if (length(v_tmp)==0) {V_sm <- rbind(V_sm,V_tmp[k,]); v_tmp <- nrow(V_sm)}
			T_tmp2 <- c(T_tmp2,v_tmp)
		  }
		  T_sm <- rbind(T_sm,T_tmp2)
		}
		upp.tri <- which(rownames(T_sm)=="T_tmp2")
	  } else { 
		low.tri = upp.tri = NULL; 
		print("We can not model the continuty in the second direction") 
	  }
	  if (flag.left) {
		left.tri <- left.tri[duplicated(left.tri)]; 
		for (j in 1: length(left.tri)) {
		  V_tmp <- V_sm[T_sm[left.tri[j],],]; V_tmp[,1] <- V_tmp[,1] + 360
		  T_tmp1 <- NULL;
		  for (k in 1:3) {
			v_tmp <- which(apply(t(V_sm)==V_tmp[k,],2,sum)==3)
			if (length(v_tmp)==0) {V_sm <- rbind(V_sm,V_tmp[k,]); v_tmp <- nrow(V_sm)}
			T_tmp1 <- c(T_tmp1,v_tmp)
		  }
		  T_sm <- rbind(T_sm,T_tmp1)
		}
	    right.tri <- which(rownames(T_sm)=="T_tmp1")
	  } else { 
		left.tri = right.tri = NULL; 
		print("We can not model the continuty in the first direction") 
	  }
	}

 	result=tdata(V_sm, T_sm)
 	E_sm=result[[1]];TE_sm=result[[2]];TV_sm=result[[3]];EV_sm=result[[4]]
 	int_sm=which(apply(TE_sm,2,sum)>1)
 	n_sm=length(int_sm)
 	N_sm=0
 	Neq_sm=0
 	result=crarrays(d_sm,r_sm)
 	I1_sm=result[[1]]
 	I2_sm=result[[2]]
  	I_sm=list();J_sm=list();K_sm=list()
 	for(j in 0:r_sm){
 		N_sm=N_sm+((j+1)*(j+2)/2+1)*(d_sm+1-j)
 		Neq_sm=Neq_sm+d_sm+1-j
 		result=indices(j)
 		LI_sm=result[,1]
 		LJ_sm=result[,2]
 		LK_sm=result[,3]
 		I_sm[[j+1]]=LI_sm
 		J_sm[[j+1]]=LJ_sm
 		K_sm[[j+1]]=LK_sm
 	}
 	nbasis_sm=(d_sm+1)*(d_sm+2)/2
 	Index1_sm=matrix(0,N_sm*n_sm,1)
 	Index2_sm=matrix(0,N_sm*n_sm,1)
 	values_sm=Index1_sm
 	A_sm=matrix(c(1,2,2,3,3,1,2,1,3,2,1,3),ncol=2,byrow=T)
 	for (j in 1:n_sm){
 		k_sm=int_sm[j]
 		v1_sm=E_sm[k_sm,1]
 		v2_sm=E_sm[k_sm,2]
 		AdjT_sm=which(TE_sm[,k_sm]!=0)
 		t1_sm=AdjT_sm[1]
 		t2_sm=AdjT_sm[2]
 		T1_sm=T_sm[AdjT_sm[1],]
 		T2_sm=T_sm[AdjT_sm[2],]
 		i1_sm=which(T1_sm==v1_sm)
 		i2_sm=which(T1_sm==v2_sm)
 		j1_sm=which(T2_sm==v1_sm)
 		j2_sm=which(T2_sm==v2_sm)
 		e1_sm = locate(matrix(c(i1_sm,i2_sm),nrow=1),A_sm)
 		e2_sm = locate(matrix(c(j1_sm,j2_sm),nrow=1),A_sm)
 		if (e1_sm >3) {
 			e1_sm=e1_sm-3
 			temp_sm=T1_sm
 			T1_sm=T2_sm
 			T2_sm=temp_sm
 			temp_sm=e1_sm
 			e1_sm=e2_sm
 			e2_sm=temp_sm
 			temp_sm=t1_sm
 			t1_sm=t2_sm
 			t2_sm=temp_sm
 		} else
 			e2_sm=e2_sm-3
 
 		v4_sm=T2_sm[T2_sm%in%c(v1_sm,v2_sm)==FALSE]
 		V4_sm=V_sm[v4_sm,]
 		result=xy2bary(V_sm[T1_sm,],V4_sm)
 		lam1_sm=result[,1]
 		lam2_sm=result[,2]
 		lam3_sm=result[,3]
 		lambda_sm=c(lam1_sm,lam2_sm,lam3_sm)
 		
 		if(e1_sm==2){
 			temp_sm=lambda_sm[1]
 			lambda_sm[1]=lambda_sm[2]
			lambda_sm[2]=lambda_sm[3]
			lambda_sm[3]=temp_sm
 		}

 		if(e1_sm==3){
 			temp_sm=lambda_sm[2]
 			lambda_sm[2]=lambda_sm[3]
 			lambda_sm[3]=temp_sm
 		}
 		VarCT_sm=0
 		EqCt_sm=0
 		for (k in 0:r_sm){
 			lam_sm=factorial(k)/(gamma(I_sm[[k+1]]+1)*gamma(J_sm[[k+1]]+1)*gamma(K_sm[[k+1]]+1))*lambda_sm[1]^I_sm[[k+1]]*lambda_sm[2]^J_sm[[k+1]]*lambda_sm[3]^K_sm[[k+1]]
 			T1mat_sm=as.matrix(I1_sm[[k+1]][[e1_sm]])
 			T2vector_sm=I2_sm[[k+1]][[e2_sm]]
 			numeq_sm=nrow(T1mat_sm)
 			numVar_sm=(ncol(T1mat_sm)+1)*numeq_sm
 			T1Values_sm=matrix(1,nrow(T1mat_sm),ncol(T1mat_sm))%*%diag(lam_sm)
 			eqnums_sm=(j-1)*(Neq_sm)+EqCt_sm+diag(1:numeq_sm)%*%matrix(1,numeq_sm,ncol(T1mat_sm)+1)
 			Index1_sm[((j-1)*N_sm+VarCT_sm+1):((j-1)*N_sm+VarCT_sm+numVar_sm)]=eqnums_sm
 			Index2_sm[((j-1)*N_sm+VarCT_sm+1):((j-1)*N_sm+VarCT_sm+numVar_sm)]=c((t1_sm-1)*nbasis_sm+T1mat_sm,(t2_sm-1)*nbasis_sm+T2vector_sm)
 			values_sm[((j-1)*N_sm+VarCT_sm+1):((j-1)*N_sm+VarCT_sm+numVar_sm)]=c(T1Values_sm,-1*rep(1,length(T2vector_sm)))
 			VarCT_sm=VarCT_sm+numVar_sm
 			EqCt_sm=EqCt_sm+numeq_sm
 		}
 	}
 	H_sm=matrix(0,n_sm*Neq_sm,nrow(T_sm)*nbasis_sm)
 	H_sm[cbind(Index1_sm,Index2_sm)]=values_sm
	
	#Fix the smoothness matrix for angular data
	if (type=="angular") {
		A <- diag((nrow(T_sm)-length(c(upp.tri,right.tri)))*nbasis_sm); A <- rbind(A,matrix(0,nrow=length(c(upp.tri,right.tri))*nbasis_sm,ncol=ncol(A)))
		remove.ind <- which(apply(abs(H_sm)[,(1:ncol(A))],1,sum)==0); if (length(remove.ind)) H_sm <- H_sm[-remove.ind,]
		indx11 <- rep(upp.tri,each=nbasis_sm)*nbasis_sm+rep(((1-nbasis_sm):0),length(upp.tri))
		indx12 <- rep(low.tri,each=nbasis_sm)*nbasis_sm+rep(((1-nbasis_sm):0),length(low.tri))
		indx21 <- rep(right.tri,each=nbasis_sm)*nbasis_sm+rep(((1-nbasis_sm):0),length(right.tri))
		indx22 <- rep(left.tri,each=nbasis_sm)*nbasis_sm+rep(((1-nbasis_sm):0),length(left.tri))
		if (length(indx11)!=1) diag(A[indx11,indx12]) <- 1 else A[indx11,indx12] <- 1
		if (length(indx21)!=1) diag(A[indx21,indx22]) <- 1 else A[indx21,indx22] <- 1
		H_sm <- H_sm%*%A
	}

	if (intercept==TRUE) {
#		ntri_sm <- nrow(mesh_sm$graph$tv);
#		H_sm <- rbind(H_sm,t(matrix(c(rep(1,nbasis_sm),rep(c(rep(0,ntri_sm*nbasis_sm),rep(1,nbasis_sm)),(ntri_sm-1))),nr=(ntri_sm*nbasis_sm))))
		H_sm <- rbind(1,H_sm)
	}
	
 	return(H_sm)
} 

################## Used in smoothness. Return the vertex, edge and etc given triangulation. #########
tdata=function(v_t,t_t){
 	ntri_t=nrow(t_t)
 	Edges_t=matrix(rep(integer(0),4),ncol=2)
 	TE_t=rep(integer(0),ntri_t)
 	numEdges_t=0
 	for (j in 1:ntri_t){
 		Tj_t=t_t[j,]
 		for (k in 1:3){
 			edge_t=c(min(Tj_t[k],Tj_t[k%%3+1]),max(Tj_t[k],Tj_t[k%%3 + 1]))
 			if (nrow(Edges_t)>0 ) {edgenum_t=which((edge_t[1]==Edges_t[,1])*(edge_t[2]==Edges_t[,2])==1)} else
 			{edgenum_t=integer(0)}
 			
 			if (length(edgenum_t)==0) {
 				Edges_t=rbind(Edges_t,edge_t)
 				numEdges_t=numEdges_t+1
 				TE_t=cbind(TE_t, rep(0,ntri_t))
 				edgenum_t=numEdges_t
 				}
 			TE_t[j,edgenum_t]=1
 			}
 		}
 	numV_t=dim(v_t)[1]
 	TV_t=matrix(0,ntri_t,numV_t)
 	for (j in 1:ntri_t){
 		TV_t[j,t_t[j,]]=1
 		}
 	EV_t=matrix(0,numEdges_t,numV_t)
 	for (j in 1:numEdges_t){
 		EV_t[j,Edges_t[j,]]=1
 		}
 	return(list(Edges_t,TE_t,TV_t,EV_t))
}

################# Used in smoothness ######
crarrays=function(d_cr,r_cr){
 	I1_cr=list();I2_cr=list()
 	for (j in  0:r_cr){	
 		result=crcellarrays(d_cr,j)
       I1_cr[[j+1]]=result[[1]]
       I2_cr[[j+1]]=result[[2]] 	
 		}
   return(list(I1_cr,I2_cr))
}
 
################# Used in smoothness ######
crcellarrays=function(d_cs,r_cs){
 	I1_cs=list();I2_cs=list();
 	result=cr_indices(d_cs,r_cs)
 	J1_cs=result[[1]];J2_cs=result[[2]]
 	D1_cs=rep(0,d_cs+1)
 	D2_cs=rep(0,d_cs+1)
 	s1_cs=d_cs+1
 	s2_cs=1
 	for (j in 1:(d_cs+1)){
 		D1_cs[j]=s1_cs
 		s1_cs=s1_cs+d_cs+1-j
 		D2_cs[j]=s2_cs
 		s2_cs=s2_cs+d_cs+2-j
 	}
 	I2_cs[[1]]=flip(J2_cs)
 	Temp_cs=D1_cs-r_cs	
 	I2_cs[[2]]=flip(Temp_cs[1:(d_cs+1-r_cs)])
 	Temp_cs=D2_cs+r_cs
 	I2_cs[[3]]=as.matrix(Temp_cs[1:(d_cs+1-r_cs)])
 	I1_cs[[1]]=J1_cs
 	Temp_cs=matrix(0,nrow(J1_cs),ncol(J1_cs))
 	for (j in 0:r_cs){Temp_cs[,j+1]=D1_cs[(j+1):(d_cs+1-r_cs+j)]}
 	loc_cs=r_cs+2
 	back_cs=r_cs+1
 	if(r_cs>0){
 	  for (j in 1:r_cs){
 		for(k in 0:(r_cs-j)){
 			Temp_cs[,loc_cs]=Temp_cs[,loc_cs-back_cs]-1
 			loc_cs=loc_cs+1
 			}
 		back_cs=back_cs-1
 	  }
 	}
 	I1_cs[[2]]=Temp_cs
 	Temp_cs=matrix(0,nrow(J1_cs),ncol(J1_cs))
 	for(j in 0:r_cs){Temp_cs[,j+1]=D2_cs[(j+1):(d_cs+1-r_cs+j)]}
 	loc_cs=r_cs+2
 	back_cs=r_cs+1
 	if(r_cs>0){
 	  for (j in 1:r_cs){
 		for(k in 0:(r_cs-j)){
 			Temp_cs[,loc_cs]=Temp_cs[,loc_cs-back_cs]+1
 			loc_cs=loc_cs+1
 		}
 		back_cs=back_cs-1
 	  }
	}
	I1_cs[[3]]=flip(Temp_cs)
	return(list(I1_cs,I2_cs))
}
 
################# Used in smoothness ######
cr_indices=function(d_ci,r_ci){
 	I1_ci=integer(0)
 	start_ci=1
 	D_ci=d_ci+1
 	for (j in 0:r_ci){
 		for (k in 0:(r_ci-j)){
 			new_col_ci=(start_ci+k):(start_ci+k+d_ci-r_ci)
 			I1_ci=cbind(I1_ci,new_col_ci)
 		}
 		start_ci=start_ci+D_ci
 		D_ci=D_ci-1
 	}
    I2_ci=(-r_ci*r_ci/2+r_ci*(d_ci+3/2)+1):(-r_ci*r_ci/2+r_ci*(d_ci+1/2)+d_ci+1)
    return(list(as.matrix(I1_ci),as.matrix(I2_ci)))
}

################# Used in smoothness ######
flip=function(matrix_fl){
 	matrix_fl=as.matrix(matrix_fl)
 	n_fl=nrow(matrix_fl)
 	m_fl=matrix(0,n_fl,n_fl)
 	for(j in 1:n_fl){m_fl[j,n_fl-j+1]=1}
 	return(m_fl%*%matrix_fl)
}
#######################################################################     END OF SMOOTHNESS        ###################################################################



#######################################################################         PENALTY MATRIX        ###################################################################

################## Numerical integration to obtain the penalty matrix: by Gaussian quadrature for general triangular elements.
int_pen=function(mesh.vert,d=3) {
	mesh.int <<- list(); zz <- 0; N=d;
	xw <- TriGaussPoints(N); 
	NP=length(xw[,1]);
	if (is.matrix(mesh.vert)) {mesh.int$loc <<- mesh.vert[,1:2]; mesh.int$graph$tv <<- matrix(1:3,nrow=1)} else {mesh.int <<- mesh.vert}
	for (i in 1:nrow(mesh.int$graph$tv)) {
		z = 0.0;
		xy <- mesh.int$loc[mesh.int$graph$tv[i,],][,1:2];
		A=abs(xy[1,1]*(xy[2,2]-xy[3,2])+xy[2,1]*(xy[3,2]-xy[1,2])+xy[3,1]*(xy[1,2]-xy[2,2]))/2;
		for (j in 1:NP) {
			x = xy[1,1]*(1-xw[j,1]-xw[j,2])+xy[2,1]*xw[j,1]+xy[3,1]*xw[j,2];
			y = xy[1,2]*(1-xw[j,1]-xw[j,2])+xy[2,2]*xw[j,1]+xy[3,2]*xw[j,2];
				b11 <- dirder2(cbind(x,y),i,d,dir1=c(0,1),dir2=c(0,1));
				b12 <- dirder2(cbind(x,y),i,d,dir1=c(0,1),dir2=c(1,0));
				b22 <- dirder2(cbind(x,y),i,d,dir1=c(1,0),dir2=c(1,0));
			f.x_y <- (t(b11)%*%b11+2*t(b12)%*%b12+t(b22)%*%b22);
			z = z + f.x_y * xw[j,3]
		}
		zz = zz + A*z; print(c(i,j));
	}
	return(zz)
}

############ Gaussian quadrature Triangles ########
TriGaussPoints=function(n) {
   if (n == 1) {
	xw=matrix(c(0.33333333333333, 0.33333333333333, 1.00000000000000),nrow=1);
   } else if (n == 2) {
	xw=matrix(c(0.16666666666667, 0.16666666666667, 0.33333333333333,
				0.16666666666667, 0.66666666666667, 0.33333333333333,
				0.66666666666667, 0.16666666666667, 0.33333333333333),ncol=3,byrow=T)
   } else if (n == 3) {
	xw=matrix(c(0.33333333333333, 0.33333333333333, -0.56250000000000,
				0.20000000000000, 0.20000000000000, 0.52083333333333,
				0.20000000000000, 0.60000000000000, 0.52083333333333,
				0.60000000000000, 0.20000000000000, 0.52083333333333),ncol=3,byrow=T)
   } else if (n == 4) {
	xw=matrix(c(0.44594849091597, 0.44594849091597, 0.22338158967801,
				0.44594849091597, 0.10810301816807, 0.22338158967801,
				0.10810301816807, 0.44594849091597, 0.22338158967801,
				0.09157621350977, 0.09157621350977, 0.10995174365532,
				0.09157621350977, 0.81684757298046, 0.10995174365532,
				0.81684757298046, 0.09157621350977, 0.10995174365532),ncol=3,byrow=T)
   } else if (n == 5) {
	xw=matrix(c(0.33333333333333, 0.33333333333333, 0.22500000000000,
				0.47014206410511, 0.47014206410511, 0.13239415278851,
				0.47014206410511, 0.05971587178977, 0.13239415278851,
				0.05971587178977, 0.47014206410511, 0.13239415278851,
				0.10128650732346, 0.10128650732346, 0.12593918054483,
				0.10128650732346, 0.79742698535309, 0.12593918054483,
				0.79742698535309, 0.10128650732346, 0.12593918054483),ncol=3,byrow=T)
   } else if (n == 6) {
	xw=matrix(c(0.24928674517091, 0.24928674517091, 0.11678627572638,
				0.24928674517091, 0.50142650965818, 0.11678627572638,
				0.50142650965818, 0.24928674517091, 0.11678627572638,
				0.06308901449150, 0.06308901449150, 0.05084490637021,
				0.06308901449150, 0.87382197101700, 0.05084490637021,
				0.87382197101700, 0.06308901449150, 0.05084490637021,
				0.31035245103378, 0.63650249912140, 0.08285107561837,
				0.63650249912140, 0.05314504984482, 0.08285107561837,
				0.05314504984482, 0.31035245103378, 0.08285107561837,
				0.63650249912140, 0.31035245103378, 0.08285107561837,
				0.31035245103378, 0.05314504984482, 0.08285107561837,
				0.05314504984482, 0.63650249912140, 0.08285107561837),ncol=3,byrow=T)
   } else if (n == 7) {
	xw=matrix(c(0.33333333333333, 0.33333333333333, -0.14957004446768,
				0.26034596607904, 0.26034596607904, 0.17561525743321,
				0.26034596607904, 0.47930806784192, 0.17561525743321,
				0.47930806784192, 0.26034596607904, 0.17561525743321,
				0.06513010290222, 0.06513010290222, 0.05334723560884,
				0.06513010290222, 0.86973979419557, 0.05334723560884,
				0.86973979419557, 0.06513010290222, 0.05334723560884,
				0.31286549600487, 0.63844418856981, 0.07711376089026,
				0.63844418856981, 0.04869031542532, 0.07711376089026,
				0.04869031542532, 0.31286549600487, 0.07711376089026,
				0.63844418856981, 0.31286549600487, 0.07711376089026,
				0.31286549600487, 0.04869031542532, 0.07711376089026,
				0.04869031542532, 0.63844418856981, 0.07711376089026),ncol=3,byrow=T)
   } else if (n == 8) {
	xw=matrix(c(0.33333333333333, 0.33333333333333, 0.14431560767779,
				0.45929258829272, 0.45929258829272, 0.09509163426728,
				0.45929258829272, 0.08141482341455, 0.09509163426728,
				0.08141482341455, 0.45929258829272, 0.09509163426728,
				0.17056930775176, 0.17056930775176, 0.10321737053472,
				0.17056930775176, 0.65886138449648, 0.10321737053472,
				0.65886138449648, 0.17056930775176, 0.10321737053472,
				0.05054722831703, 0.05054722831703, 0.03245849762320,
				0.05054722831703, 0.89890554336594, 0.03245849762320,
				0.89890554336594, 0.05054722831703, 0.03245849762320,
				0.26311282963464, 0.72849239295540, 0.02723031417443,
				0.72849239295540, 0.00839477740996, 0.02723031417443,
				0.00839477740996, 0.26311282963464, 0.02723031417443,
				0.72849239295540, 0.26311282963464, 0.02723031417443,
				0.26311282963464, 0.00839477740996, 0.02723031417443,
				0.00839477740996, 0.72849239295540, 0.02723031417443),ncol=3,byrow=T)
   } else {
	print("Bad input n");
   }
   return(xw)
}

########## Second derivatives used in the penalty ###########
dirder2=function(data,i,d=3,dir1=c(1,0),dir2=c(1,0)) {
  r=2; ind=indices(d); ind.d_r=indices(d-r)
  n.basis=(d+1)*(d+2)/2; n.basis.d_r=(d-r+1)*(d-r+2)/2
  b.matrix=matrix(0,nrow(data),n.basis*nrow(mesh.int$graph$tv))
  b.d_r=basis(mesh.int,d-r,data)
	for (j in 1:n.basis){
		tdir1=as.vector(dcord(mesh.int$loc[mesh.int$graph$tv[i,],],dir1)); tdir2=as.vector(dcord(mesh.int$loc[mesh.int$graph$tv[i,],],dir2));
		tdir=c(tdir1*tdir2,tdir1[c(1,1,2)]*tdir2[c(2,3,3)]+tdir2[c(1,1,2)]*tdir1[c(2,3,3)])
		ind.6 <- rbind(ind[j,]-c(2,0,0),ind[j,]-c(0,2,0),ind[j,]-c(0,0,2),ind[j,]-c(1,1,0),ind[j,]-c(1,0,1),ind[j,]-c(0,1,1));
		ind.rm <- which(ind.6<0,arr.ind=T); ind.r <- vector()
		if (length(ind.rm)) {ind.6 <- matrix(ind.6[-ind.rm[,1],],nc=3); tdir <- tdir[-ind.rm[,1]]}
		for(k in 1:nrow(ind.6)){ind.r <- c(ind.r, which(apply(ind.d_r, 1, function(x) all(x == ind.6[k,]))))}
		b.matrix[,(i-1)*n.basis+j] = apply(t(tdir*t(matrix(b.d_r[,(i-1)*n.basis.d_r+ind.r],nr=nrow(data)))),1,sum)
	}
  return(b.matrix)
}

########## directional coordinates used in the derivatives ##########
dcord=function(triangle,v_tc){
	result1=xy2bary(triangle,v_tc)
	result2=xy2bary(triangle,c(0,0))
	return(result1-result2)
}

########## First derivatives (Not used anywhere) ###########
dirder1=function(data,mesh,d=3,dir=c(1,0)) {
  r=1; ind=indices(d); ind.d_r=indices(d-r)
  n.basis=(d+1)*(d+2)/2; n.basis.d_r=(d-r+1)*(d-r+2)/2
  b.matrix=matrix(0,nrow(data),n.basis*nrow(mesh$graph$tv))
  b.d_r=basis(mesh,d-r,data)
  for (i in 1:nrow(mesh$graph$tv)) {
    result <- xy2bary(mesh$loc[mesh.$graph$tv[i,],],data)
    indx <- apply(round(result,15),1,min)>=0
    if (sum(indx)) {
      for (j in 1:n.basis) {
		tdir=as.vector(dcord(mesh$loc[mesh$graph$tv[i,],],dir))
		ind.3 <- rbind(ind[j,]-c(1,0,0),ind[j,]-c(0,1,0),ind[j,]-c(0,0,1));
		ind.rm <- which(ind.3<0,arr.ind=T); ind.r <- vector()
		if (length(ind.rm)) {ind.3 <- matrix(ind.3[-ind.rm[,1],],nc=3); tdir <- tdir[-ind.rm[,1]]}
		for(k in 1:nrow(ind.3)){ind.r <- c(ind.r, which(apply(ind.d_r, 1, function(x) all(x == ind.3[k,]))))}
#		ind[j,];ind.d_r[ind.r,];tdir
		b.matrix[indx,(i-1)*n.basis+j] = apply(t(tdir*t(matrix(b.d_r[indx,(i-1)*n.basis.d_r+ind.r],nr=sum(indx)))),1,sum)
	  }
    }
  }
  return(b.matrix)
}



################################################# Huijun's code for calculating the penalty (Much faster, works fine just for d=3, r=1) ##############################
  ###################### energy
  ###################### evaluate energy matrix K
  energy=function(mesh_e,d_e,index_e){
  	V_e <- mesh_e$loc; T_e <- mesh_e$graph$tv;

 	ntri_e=nrow(T_e)
 	D_e=(d_e+1)*(d_e+2)/2
 	Dsq_e=D_e^2
 	Mat_e=build(d_e-2)
 	Index1_e=rep(0,ntri_e*Dsq_e)
 	Index2_e=Index1_e
 	S_e=Index1_e
 	place_e=1
 	for (k in 1:ntri_e){
 		LocK_e=locEng(V_e[T_e[k,1],],V_e[T_e[k,2],],V_e[T_e[k,3],],Mat_e,d_e,index_e)
 		result=which(LocK_e!=0,arr.in=TRUE)
 		i_e=result[,1];j_e=result[,2];s_e=LocK_e[cbind(i_e,j_e)]
  		L_e=length(i_e)
 		Index1_e[place_e:(place_e+L_e-1)]=(k-1)*D_e+i_e
 		Index2_e[place_e:(place_e+L_e-1)]=(k-1)*D_e+j_e
 		S_e[place_e:(place_e+L_e-1)]=s_e
 		place_e=place_e+L_e
 		}
		K_e=matrix(0,ntri_e*D_e,ntri_e*D_e)
 		K_e[cbind(Index1_e[1:(place_e-1)],Index2_e[1:(place_e-1)])]=S_e[1:(place_e-1)]
 		return(K_e)
 	}
 	
  ######################## locEng
  ######################## calculate energy on each triangle
  locEng=function(v1_lo,v2_lo,v3_lo,Mat_lo,d_lo,index_lo){
	D_lo=(d_lo+1)*(d_lo+2)/2
	Id_lo=diag(rep(1,D_lo))
	vx_lo=c(1,0)
	vy_lo=c(0,1)
	result=tcord(v1_lo,v2_lo,v3_lo,vx_lo)
	lam1x_lo=result[,1];lam2x_lo=result[,2];lam3x_lo=result[,3]
	result=tcord(v1_lo,v2_lo,v3_lo,vy_lo)
	lam1y_lo=result[,1];lam2y_lo=result[,2];lam3y_lo=result[,3]
	
	Dx_lo=dirder(Id_lo,lam1x_lo,lam2x_lo,lam3x_lo)
	Dxx_lo=dirder(Dx_lo,lam1x_lo,lam2x_lo,lam3x_lo)
	Dxy_lo=dirder(Dx_lo,lam1y_lo,lam2y_lo,lam3y_lo)
	Dy_lo=dirder(Id_lo,lam1y_lo,lam2y_lo,lam3y_lo)
	Dyy_lo=dirder(Dy_lo,lam1y_lo,lam2y_lo,lam3y_lo)

	if(index_lo==1){
	K3_lo=abs(triarea(v1_lo,v2_lo,v3_lo))*(t(Dxx_lo)%*%Mat_lo%*%Dxx_lo+2*t(Dxy_lo)%*%Mat_lo%*%Dxy_lo+t(Dyy_lo)%*% Mat_lo%*%Dyy_lo)
	}
	if(index_lo==2){
	K3_lo=abs(triarea(v1_lo,v2_lo,v3_lo))*(t(Dxx_lo+Dyy_lo)%*%Mat_lo%*%(Dxx_lo+Dyy_lo))
	}
	return(K3_lo)
	}

  ################## triarea
  ################## calcuate the area of triangle given vertice
  triarea=function(v1_ta,v2_ta,v3_ta){
 	x_ta=v1_ta[1];y_ta=v1_ta[2]
 	a_ta=v2_ta[1];b_ta=v2_ta[2]
 	c_ta=v3_ta[1];d_ta=v3_ta[2]
 	A_ta=(a_ta-x_ta)*(d_ta-y_ta)-(c_ta-x_ta)*(b_ta-y_ta)
 	A_ta=A_ta/2
 	return(A_ta)
 	}
 
  ######################## tcord
  ######################## return the barycentric coord of direction
  tcord=function(v1_tc,v2_tc,v3_tc,v_tc){
	result1=xy2bary(rbind(v1_tc,v2_tc,v3_tc),v_tc)
	result2=xy2bary(rbind(v1_tc,v2_tc,v3_tc),c(0,0))
	return(result1-result2)
	}
 	
  ######################## build
  ######################## evaluate the matrix for inner product
  build=function(d_bu){
	result=indices(d_bu)
	I_bu=result[,1];J_bu=result[,2];K_bu=result[,3]
	m_bu=(d_bu+1)*(d_bu+2)/2
	Mat_bu=matrix(0,m_bu,m_bu)
	for (j in 1:m_bu){
		for(k in 1:m_bu){
			Mat_bu[k,j]=choose(I_bu[j]+I_bu[k],I_bu[j])*choose(J_bu[j]+J_bu[k],J_bu[j])*choose(K_bu[j]+K_bu[k],K_bu[j])
			}
		}
	Mat_bu=Mat_bu/(choose(d_bu*2,d_bu)*choose(2*d_bu+2,2))
	return(Mat_bu)
	}
	
  ####################### dirder
  ####################### 
  dirder=function(Bc_di,lam1_di,lam2_di,lam3_di){
	m_di=dim(Bc_di)[1]
	d_di=degree(m_di)
	DerBc_di=d_di*de_cast_step(lam1_di,lam2_di,lam3_di,Bc_di)
	return(DerBc_di)
	}
#######################################################################      END OF PENALTY MATRIX    ###################################################################



#######################################################################         REST OF ...           ###################################################################
	
################# The Simulation in the paper ###########
simulate <- function(nclass=20,nobs=50,x.range=c(-180,180),y.range=c(-180,180),n.grid=45,setup=1) {
  angs <- list(); dens <- list(); dens.o <- list();
  if (setup==1) {p <- c(0.04,0.2,0.96)} else {p <- c(0.04,0.96)};
  ps <- sample(rbind(1:length(p),1:nclass)[1,]);
  if (setup==1) {
	theta.p <- as.circular(c(95,120),units="degree"); tau.p <- as.circular(c(50,-165),units="degree"); theta.r=c(0.995,0.99); tau.r=c(0.975,0.975)
  } else if (setup==2) {
	theta.p <- as.circular(c(95,135),units="degree"); tau.p <- as.circular(c(95,-135),units="degree"); theta.r=c(0.995,0.99); tau.r=c(0.975,0.975)
  } else{
	mu <- Sigma <- list(); mu[[1]] <- c(95,95); mu[[2]] <- c(135,-135); Sigma[[1]] <- diag(c(40,160)); Sigma[[2]] <- matrix(c(70,-20,-20,140),nr=2);
  }
  if (setup!=3) {
	xs <- as.circular(seq(x.range[1],x.range[2],length.out=n.grid),units="degree");
	ys <- as.circular(seq(y.range[1],y.range[2],length.out=n.grid),units="degree")
	for (i in 1:2) dens.o[[i]] <- outer(dwrappednormal(xs, mu=theta.p[i], rho=theta.r[i]), dwrappednormal(ys, mu=tau.p[i], rho=tau.r[i]), "*")
  } else {
	xs <- seq(x.range[1],x.range[2],length.out=n.grid); ys <- seq(y.range[1],y.range[2],length.out=n.grid)
	for (i in 1:2) dens.o[[i]] <- matrix(dmvnorm(cbind(rep(xs,length(ys)),rep(ys,each=length(xs))),mean=mu[[i]],sigma=Sigma[[i]]),nr=n.grid)
  }
	for (i in 1:length(p)) {
	  dens[[i]] <- list(); dens[[i]]$x <- as.vector(xs); dens[[i]]$y <- as.vector(ys); 
	  dens[[i]]$z <- array(0, dim=c(length(xs),length(ys))); dens[[i]]$z <- p[i]*dens.o[[1]]+(1-p[i])*dens.o[[2]]
	  if (setup!=3) dens[[i]]$z <- dens[[i]]$z*(pi/180)^2
	}
	for (i in 1:nclass) {
	  if (setup!=3) {
		theta <- c(rwrappednormal(p[ps[i]]*nobs, mu=theta.p[1], rho=theta.r[1]),rwrappednormal((1-p[ps[i]])*nobs, mu=theta.p[2], rho=theta.r[2]))
		tau <- c(rwrappednormal(p[ps[i]]*nobs, mu=tau.p[1], rho=tau.r[1]),rwrappednormal((1-p[ps[i]])*nobs, mu=tau.p[2], rho=tau.r[2]))
		angs[[i]] <- cbind(theta,tau); angs[[i]][angs[[i]]>180] <- angs[[i]][angs[[i]]>180]-360
	  } else {
		angs[[i]] <- rbind(rmvnorm(round(p[ps[i]]*nobs), mean=mu[[1]], sigma=Sigma[[1]]),rmvnorm(round((1-p[ps[i]])*nobs), mean=mu[[2]], sigma=Sigma[[2]]))
	  }
	}
  return(list(angs=angs,cols=ps,dens=dens))
}

################## Kernel Density Estimates (Angular) ###########
kde2d <- function (x, y, h, n = 25, lims = c(range(x), range(y)), polar=FALSE) {
    nx <- length(x)
    if (length(y) != nx) stop("data vectors must be the same length")
    if (any(!is.finite(x)) || any(!is.finite(y))) stop("missing or infinite values in the data are not allowed")
    if (any(!is.finite(lims))) stop("only finite values are allowed in 'lims'")
    n <- rep(n, length.out = 2L)
    gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
    gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
    h <- if (is.null(h)) c(bandwidth.nrd(x), bandwidth.nrd(y))
    else rep(h, length.out = 2L)
    h <- h/4;
    ax <- outer(gx, x, "-")
    ay <- outer(gy, y, "-")
    if (polar) {
        ax[ax>180] <- 360-ax[ax>180]
        ax[ax<(-180)] <- 360+ax[ax<(-180)]
        ay[ay>180] <- 360-ay[ay>180]
        ay[ay<(-180)] <- 360+ay[ay<(-180)]
    }
    ax <- ax/h[1L]
    ay <- ay/h[2L]
    z <- tcrossprod(matrix(dt(ax,20), , nx), matrix(dt(ay,20), , nx))/(nx * h[1L] * h[2L])
    list(x = gx, y = gy, z = z)
}

################## Set of KDEs ###########
densities <- function(angs,x1.range=c(-180,180),x2.range=c(-180,180),n.grid=45,type="angular",sqrt.t=FALSE,h=NULL,res=NULL) {
	Xs <- list(); len=length(angs); if (type=="angular") {polar=TRUE} else {polar=FALSE}
	for (i in 1:len) {
	  Xs[[i]] <- kde2d(angs[[i]][,1],angs[[i]][,2],n=n.grid,lims=c(x1.range,x2.range),polar=polar,h=h)
#	  if (!is.null(res)) Xs[[i]]$z <- (Xs[[i]]$z + min(res$dens[[cols[i]]]$z)); Xs[[i]]$z <- Xs[[i]]$z/sum(Xs[[i]]$z*temp)
	  Xs[[i]]$z[which((Xs[[i]]$z)==0)]=5e-324
	  if (sqrt.t) Xs[[i]]$z <- sqrt(Xs[[i]]$z); 
	}
	return(Xs)
}

dens.bnd <- function(angs,x1.range=c(-180,180),x2.range=c(-180,180),n.grid=45,sqrt.t=FALSE,methods=c("Hpi","Hbcv","Hscv","Hlscv"),res=NULL) {
	Xs <- list(); len=length(angs);
	for (j in 1:length(methods)) {
	 Xs[[j]] <- list()
	 for (i in 1:len) {
	  if (methods[j]=="Hpi") bnd <- Hpi(angs[[i]])
	  if (methods[j]=="Hbcv") bnd <- Hbcv(angs[[i]])
	  if (methods[j]=="Hscv") bnd <- Hscv(angs[[i]])
	  if (methods[j]=="Hlscv") bnd <- Hlscv(angs[[i]])
	  tmp <- kde(angs[[i]],bnd,gridsize=n.grid,xmin=c(x1.range,x2.range)[c(1,3)],xmax=c(x1.range,x2.range)[c(2,4)])
#	  if (!is.null(res)) tmp$estimate <- (tmp$estimate- min(tmp$estimate)) * diff(range(res$dens[[cols[i]]]$z)) / diff(range(tmp$estimate)) + min(res$dens[[cols[i]]]$z)
#	  if (!is.null(res)) tmp$estimate <- (tmp$estimate + min(res$dens[[cols[i]]]$z)); tmp$estimate <- tmp$estimate/sum(tmp$estimate*temp)
	  Xs[[j]][[i]] <- list(x=tmp$eval.points[[1]], y=tmp$eval.points[[2]], z=tmp$estimate)
	  Xs[[j]][[i]]$z[which((Xs[[j]][[i]]$z)<1e-35)]=1e-35
	  if (sqrt.t) Xs[[j]][[i]]$z <- sqrt(Xs[[j]][[i]]$z); 
	 }
	}
	return(Xs)
}

######################## Clustering the raw data to different sets ########
cluster <- function(pdb.x,ang,method="lag",len=2) {
	angs <- list();
	if (method %in% c("lag","protein")) {
		indx <- 1:nrow(pdb.x);	#indx <- which(pdb,x[,5]=="0"); ## 1:alpha, 2:beta, 0:coil
		list.n <- apply(cbind(substring(pdb.x[,1],1,4),pdb.x[,3]),1,concat)
		list.e <- c((sort(match(unique(list.n),(list.n)))-1)[-1],length(list.n))
	}
	if (method=="lag") {
		for (i in 1:len) {
			angs.1 <- lags(ang[,1],list.e,i,indx);
			angs.2 <- lags(ang[,2],list.e,i,indx);
			angs[[i]] <- cbind(angs.1[,1],angs.2[,2]);
		}
	} else if (method=="protein") {
		list.s <- c(0,(list.e[-length(list.e)]))+1;
		for (i in 1:max(length(list.e),len)) {
			angs[[i]] <- cbind(ang[list.s[i]:list.e[i],1],ang[list.s[i]:list.e[i],2]);
		}
	} else {
		list.a <- c("GLN","ASN","LEU","VAL","SER","ILE","ASP","LYS","THR","GLU","PHE","ALA","MET","PRO","GLY","ARG","TYR","HIS","CYS","TRP");
		for (i in 1:max(length(list.a),len)) {
			list.ind <- which(pdb.x[,6]==list.a[i])
			angs[[i]] <- cbind(ang[list.ind,1],ang[list.ind,2]);
		}
	}
	return(angs)
}

################## Used in cluster ##############
concat <- function(v) {
  res = "";
  for (i in 1:length(v)) res = paste(res,v[i],sep="")
  res
}

################## Used in cluster ##############
lags <- function(angles,list.e,lag,indx) {
	list.s = c(0,(list.e[-length(list.e)])) + 1; lags <- NULL;
	for (i in 1:length(list.e)) {
	    temp <- list.s[i]:(list.e[i]-lag); temp <- temp[which(temp%in%indx)]
	    if (length(temp)) lags <- rbind(lags,cbind(angles[temp],angles[(temp+lag)]));
	}
	return(lags)
}

################## Animation plots #######
ani.plots <- function(web=TRUE){
    ani.options(outdir = getwd());
    if (web) {
		saveGIF(plots(x.lab,y.lab,"imageR"), movie.name = "imageR.gif", interval = .25, nmax = 50)
		saveGIF(plots(x.lab,y.lab,"imageS"), movie.name = "imageS.gif", interval = .25, nmax = 50)
		saveGIF(plots(x.lab,y.lab,"heatR"), movie.name = "heatR.gif", interval = .25, nmax = 50)
		saveGIF(plots(x.lab,y.lab,"heatS"), movie.name = "heatS.gif", interval = .25, nmax = 50)
		saveGIF(plots(x.lab,y.lab,"perspR",cols=cols), movie.name = "perspR.gif", interval = .25, nmax = 50)
		saveGIF(plots(x.lab,y.lab,"perspS",cols=cols), movie.name = "perspS.gif", interval = .25, nmax = 50)
    } else {
	  ani.options(qpdf = "C:/Dropbox/Util/qpdf/bin/qpdf.exe")
	  saveLatex({par(mar = c(3, 3, 1, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1)
	    plots(x.lab,y.lab,"imageR")}, img.name = "imageR", ani.opts = "controls,loop,width=0.95\\textwidth", 
	   interval = 0.20, nmax = 50, ani.dev = "pdf", ani.type = "pdf", ani.width = 7, ani.height = 7, 
	   documentclass = paste("\\documentclass{article}", "\\usepackage[papersize={7in,7in},margin=0.3in]{geometry}", sep = "\n"))
	  saveLatex({par(mar = c(3, 3, 1, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1)
	    plots(x.lab,y.lab,"imageS")}, img.name = "imageS", ani.opts = "controls,loop,width=0.95\\textwidth", 
	   interval = 0.20, nmax = 50, ani.dev = "pdf", ani.type = "pdf", ani.width = 7, ani.height = 7, 
	   documentclass = paste("\\documentclass{article}", "\\usepackage[papersize={7in,7in},margin=0.3in]{geometry}", sep = "\n"))
	  saveLatex({par(mar = c(3, 3, 1, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1)
	    plots(x.lab,y.lab,"heatR")}, img.name = "heatR", ani.opts = "controls,loop,width=0.95\\textwidth", 
	   interval = 0.20, nmax = 50, ani.dev = "pdf", ani.type = "pdf", ani.width = 7, ani.height = 7, 
	   documentclass = paste("\\documentclass{article}", "\\usepackage[papersize={7in,7in},margin=0.3in]{geometry}", sep = "\n"))
	  saveLatex({par(mar = c(3, 3, 1, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1)
	    plots(x.lab,y.lab,"heatS")}, img.name = "heatS", ani.opts = "controls,loop,width=0.95\\textwidth", 
	   interval = 0.20, nmax = 50, ani.dev = "pdf", ani.type = "pdf", ani.width = 7, ani.height = 7, 
	   documentclass = paste("\\documentclass{article}", "\\usepackage[papersize={7in,7in},margin=0.3in]{geometry}", sep = "\n"))
	  saveLatex({par(mar = c(3, 3, 1, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1)
	    plots(x.lab,y.lab,"perspR",cols=cols)}, img.name = "perspR", ani.opts = "controls,loop,width=0.95\\textwidth", 
	   interval = 0.20, nmax = 50, ani.dev = "pdf", ani.type = "pdf", ani.width = 7, ani.height = 7, 
	   documentclass = paste("\\documentclass{article}", "\\usepackage[papersize={7in,7in},margin=0.3in]{geometry}", sep = "\n"))
	  saveLatex({par(mar = c(3, 3, 1, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1)
	    plots(x.lab,y.lab,"perspS",cols=cols)}, img.name = "perspS", ani.opts = "controls,loop,width=0.95\\textwidth", 
	   interval = 0.20, nmax = 50, ani.dev = "pdf", ani.type = "pdf", ani.width = 7, ani.height = 7, 
	   documentclass = paste("\\documentclass{article}", "\\usepackage[papersize={7in,7in},margin=0.3in]{geometry}", sep = "\n"))
    }
}

################# Used in ani.plots ######
plots <- function(x1="theta",x2="tau",type="perspR",aspect1=FALSE,cols=cols) {
   par(mgp=c(1.5,.5,0),mar=c(3,3,1.5,1.75),cex=1.5);
   len <- length(Xs); y.lab <- "density"; set.seed(1); j <- 0;
   if (aspect1) {X.lim=Y.lim=range(c(Xs[[1]]$x,Xs[[1]]$y))} else {X.lim=range(Xs[[1]]$x); Y.lim=range(Xs[[1]]$y)}
   if (substr(type,1,nchar(type)-1)!="persp") {is <- (1:len); js <- sample(1:len)} else {is <- sample(1:len); js <- (1:len)}
   for (i in is) {
	XS <- Xs[[i]]; j <- j + 1;
	if (substr(type,nchar(type),nchar(type))=="S") {XS$z <- exp(matrix(init$BC%*%(ans$Theta%*%t(ans$A))[,i],nr=n.grid))}
	if (substr(type,1,nchar(type)-1)=="image") image(XS,col=rainbow(100, start=.5, end=.1),main=paste("Protein",js[j]),xlab=x1,ylab=x2,xlim=X.lim,ylim=Y.lim)
	if (substr(type,1,nchar(type)-1)=="heat") print(levelplot(XS$z,row.values=XS$x,column.values=XS$y,cuts=15,xlab=x1,ylab=x2,main=paste("Protein",js[j]),contour=T,colorkey=F,aspect="fill"))
	if (substr(type,1,nchar(type)-1)=="persp") {
		z.lim <- c(0,max(XS$z)); z <- XS$z; nrz <- nrow(z); ncz <- ncol(z); 
		jet.colors <- colorRampPalette(c(8, cols[i]+1) )
		nbcol <- 12; color <- jet.colors(nbcol);
		zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
		facetcol <- cut(zfacet, nbcol)
		persp(XS,theta=50,zlim=z.lim,main=paste("Protein",js[j]),xlab=x1,ylab=x2,xlim=X.lim,ylim=Y.lim, zlab=y.lab, col=color[facetcol])
	}
   }
}
