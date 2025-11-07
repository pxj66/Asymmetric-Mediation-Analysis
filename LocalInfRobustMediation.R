MeanCov = function(Z){
  n = nrow(Z)
  p = ncol(Z)
  
  zbar = t(Z) %*% matrix(1/n, n, 1) 
  S = t(Z) %*% (diag(n) - matrix(1/n, n, n)) %*% Z / n 
  
  Results = list(zbar = zbar, S = S)
  return(Results)
} 


vech = function(A){
  l = 0
  p = nrow(A)
  ps = p * (p + 1)/2
  vA = matrix(0, ps, 1) 
  for (i in 1:p) {
    for (j in i:p) {
      l = l + 1
      vA[l, 1] = A[j, i]			
    }	
  }
  
  return(vA)	
}


Dp = function(p){
  p2 = p * p
  ps = p * (p + 1) / 2
  Dup = matrix(0, p2, ps)
  count=0
  for (j in 1:p){
    for (i in j:p){
      count = count + 1
      if (i==j){
        Dup[(j-1)*p+j, count]=1
      }else{
        Dup[(j-1)*p+i, count]=1
        Dup[(i-1)*p+j, count]=1
      }
    }
  }	
  
  return(Dup)	
}

getSE=function(theta,Omega,n){
  
  hdot=gethdot(theta)
  COV=hdot%*%(Omega/n)%*%t(hdot)  # delta method 
  se.v=sqrt(diag(COV)) # se.v of theta     
  
  a=theta[1]
  b=theta[2]
  SobelSE=sqrt(a^2*COV[2,2]+b^2*COV[1,1])
  
  se.v=c(se.v,SobelSE) # including Sobel SE
  
  return(se.v) 
  
}

gethdot=function(theta){
  
  p=3 
  ps=p*(p+1)/2 
  q=ps
  
  a=theta[1]
  b=theta[2]
  c=theta[3]
  vx=theta[4]
  vem=theta[5]
  vey=theta[6]
  
  sigmadot=matrix(NA,ps,q)
  sigmadot[1,]=c(0,0,0,1,0,0)
  sigmadot[2,]=c(vx,0,0,a,0,0)
  sigmadot[3,]=c(b*vx,a*vx,vx,a*b+c,0,0)
  sigmadot[4,]=c(2*a*vx,0,0,a^2,1,0)
  sigmadot[5,]=c((2*a*b+c)*vx,a^2*vx+vem,a*vx,a^2*b+a*c,b,0)
  sigmadot[6,]=c((2*b*c+2*a*b^2)*vx,(2*c*a+2*a^2*b)*vx+2*b*vem,(2*a*b+2*c)*vx,c^2+2*c*a*b+a^2*b^2,b^2,1)
  
  hdot=solve(sigmadot)
  
  return(hdot)
}

MLEst=function(S){
  ahat=S[1,2]/S[1,1]
  vx=S[1,1]
  # M on X
  Sxx=S[1:2,1:2]
  sxy=S[1:2,3]
  vem=S[2,2]-S[2,1]*S[1,2]/S[1,1]
  
  # Y on X and M
  invSxx=solve(Sxx)
  beta.v=invSxx%*%sxy # chat, bhat
  vey=S[3,3]-t(sxy)%*%invSxx%*%sxy
  thetaMLE=c(ahat,beta.v[2],beta.v[1],vx,vem,vey)
  return(thetaMLE)
  
}


SEML=function(Z,thetaMLE){
  n=nrow(Z)
  p=ncol(Z)
  ps=p*(p+1)/2
  q=ps
  zbar=MeanCov(Z)$zbar
  S=MeanCov(Z)$S
  Dup=Dp(p)
  InvS=solve(S)
  W=0.5*t(Dup)%*%(InvS%x%InvS)%*%Dup
  OmegaInf=solve(W)  # only about sigma, not mu
  
  
  # Sandwich-type Omega
  S12=matrix(0,p,ps)
  S22=matrix(0,ps,ps)
  
  for (i in 1:n){
    zi0=Z[i,]-zbar
    difi=zi0%*%t(zi0)-S 
    vdifi=vech(difi)   	 
    
    S12=S12+zi0%*%t(vdifi)	
    S22=S22+vdifi%*%t(vdifi)		
  }
  
  OmegaSW=S22/n # only about sigma, not mu
  
  SEinf=getSE(thetaMLE,OmegaInf,n)
  SEsw=getSE(thetaMLE,OmegaSW,n)
  
  Results=list(inf=SEinf,sand=SEsw)
  return(Results)
  
}


HuberTun=function(kappa,p){
  prob=1-kappa  
  chip=qchisq(prob,p)
  r=sqrt(chip)
  tau=(p*pchisq(chip,p+2)+ chip*(1-prob))/p	
  Results=list(r=r,tau=tau)
  return(Results)	
}


robEst=function(Z, r, tau, ep){
  
  p = ncol(Z)
  n = nrow(Z)
  
  # Starting values 	
  mu0 = MeanCov(Z)$zbar 
  Sigma0 = MeanCov(Z)$S
  Sigin = solve(Sigma0)
  
  diverg = 0 # convergence flag
  
  for (k in 1:50) {
    sumu1 = 0
    mu = matrix(0,p,1)
    Sigma = matrix(0,p,p)
    d = rep(NA,n) 
    u1 = rep(NA,n)
    u2 = rep(NA,n)
    
    for (i in 1:n) {
      zi = Z[i,] 
      zi0 = zi - mu0 
      di2 = t(zi0)%*% Sigin %*% zi0
      di = as.numeric(sqrt(di2))
      d[i] = di
      
      #get u1i,u2i
      if (di<=r) {
        u1i = 1.0
        u2i = 1.0 / tau 
      }else { 
        u1i = r / di 
        u2i = u1i^2 / tau 
      }
      u1[i] = u1i
      u2[i]=u2i
      
      sumu1 = sumu1 + u1i
      mu = mu+u1i*zi
      Sigma = Sigma + u2i * zi0 %*% t(zi0)
      
    } # end of loop i 
    
    mu1 = mu/sumu1
    Sigma1 = Sigma/n 
    Sigdif = Sigma1-Sigma0 
    dt = sum(Sigdif^2)
    
    mu0 = mu1 
    Sigma0 = Sigma1
    Sigin = solve(Sigma0)
    
    if (dt<ep) {break}
    
  } # end of loop k
  
  
  if (k == 50) {
    diverg = 1
    mu0 = rep(0, p)
    sigma0 = matrix(NA, p ,p)
    
  }
  
  theta=MLEst(Sigma0)
  
  Results=list(mu=mu0,Sigma=Sigma0,theta=theta,d=d,u1=u1,u2=u2,diverg=diverg)   
  return(Results)              
}


SErob=function(Z,mu,Sigma,theta,d,r,tau){
  n=nrow(Z)
  p=ncol(Z)
  ps=p*(p+1)/2
  q=6
  Dup=Dp(p)
  DupPlus=solve(t(Dup)%*%Dup)%*%t(Dup)
  
  InvSigma=solve(Sigma)
  sigma=vech(Sigma)
  W=0.5*t(Dup)%*%(InvSigma%x%InvSigma)%*%Dup
  
  Zr=matrix(NA,n,p) # robustly transformed data  
  A=matrix(0,p+q,p+q) 
  B=matrix(0,p+q,p+q)
  sumg=rep(0,p+q)
  
  for (i in 1:n) {
    zi=Z[i,]
    zi0=zi-mu
    di=d[i]
    
    if (di<=r) {
      u1i=1.0
      u2i=1.0/tau
      du1i=0
      du2i=0
    }else {
      u1i=r/di
      u2i=u1i^2/tau
      du1i=-r/di^2
      du2i=-2*r^2/tau/di^3
    }
    
    #robust transformed data 
    Zr[i,]=sqrt(u2i)*t(zi0)
    
    ####	gi 
    
    g1i=u1i*zi0	# defined in (24)
    vTi=vech(zi0%*%t(zi0))
    g2i=u2i*vTi-sigma	# defined in (25)
    gi=rbind(g1i,g2i)
    sumg=gi+sumg
    
    B=B+gi%*%t(gi)
    
    ####	gdoti
    
    #	derivatives of di with respect to mu and sigma
    ddmu=-1/di*t(zi0)%*%InvSigma
    ddsigma=-t(vTi)%*%W/di
    
    #	
    dg1imu=-u1i*diag(p)+du1i*zi0%*%ddmu
    dg1isigma=du1i*zi0%*%ddsigma
    dg2imu=-u2i*DupPlus%*%(zi0%x%diag(p)+diag(p)%x%zi0)+du2i*vTi%*%ddmu
    dg2isigma=du2i*vTi%*%ddsigma-diag(q)
    
    dgi=rbind(cbind(dg1imu,dg1isigma),cbind(dg2imu,dg2isigma))
    A=A+dgi
  } # end of loop n
  
  A=-1*A/n
  B=B/n       
  invA=solve(A)
  OmegaSW=invA%*%B%*%t(invA)
  OmegaSW=OmegaSW[(p+1):(p+q),(p+1):(p+q)]
  
  
  SEsw=getSE(theta,OmegaSW,n)
  SEinf=SEML(Zr,theta)$inf
  
  Results=list(inf=SEinf,sand=SEsw,Zr=Zr)  	      
  return(Results)
  
}


# source("DGPs.R")
# Data = DGP("m", "laplace", 200)
# Z <- as.matrix(Data)
# 
# MeanCov(Z)
# TunRes <- HuberTun(0.05, p = ncol(Z))
# EstRes <- robEst(Z, r = TunRes$r, tau = TunRes$tau, ep = 1e-5)
# SERes <- SErob(Z = Z, mu = EstRes$mu, Sigma = EstRes$Sigma, theta = EstRes$theta, d = EstRes$d,
#       r = TunRes$r, tau = TunRes$tau)
# 
# sqrt(EstRes$theta[1]^2 * SERes$inf[2]^2 + EstRes$theta[2]^2 * SERes$inf[1]^2)
# SERes$inf

################################################################################
LocalInf=function(Z){
  p=ncol(Z) 
  q=p*(p+1)/2
  n=nrow(Z) 
  zbar=MeanCov(Z)$zbar 	  # sample mean 
  S=MeanCov(Z)$S 		     # S: sample covariance
  invS=solve(S)
  Dup=Dp(p)
  
  What=0.5*t(Dup)%*%(invS%x%invS)%*%Dup     # defined in (6)
  T=solve(What)/n				    # defined in (15)
  
  
  Upsilon=matrix(NA,q,n)
  
  for (i in 1:n){
    zi0=Z[i,]-zbar
    Upsilon[,i]=What%*%vech(zi0%*%t(zi0)-S)  # defined in (16)
  } # end of loop i 
  
  
  # Diagnostic measures
  
  temp=t(Upsilon)%*%T%*%Upsilon
  temp2=temp%*%temp
  deno=sqrt(sum(diag((temp2)))) 
  B.v=diag(temp/deno) 		  # Diagnostica measures as proposed by Poon&Poon(1999)     
  
  # The largest conformal normal curvature 
  
  # evv=eigen(temp)
  # values=evv$values
  # vectors=evv$vectors
  # 
  # if (values[1]>abs(values[n])){
  #   lmax=vectors[,1]
  # }else{
  #   lmax=vectors[,n]
  # }
  # 
  # CNCmax=t(lmax)%*%temp%*%lmax/deno
  
  
  # Results=list(B.v=B.v,CNCmax=CNCmax)
  Results=list(B.v=B.v)
  return(Results)
  
} # end of function


# LIRes <- LocalInf(Z)
# Z[-c(order(LIRes$B.v, decreasing = T)[1: (0.05 * n)]), ]








