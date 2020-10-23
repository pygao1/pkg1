########################
# Functions
########################
ubest = function(Y,swmat){
# compute permutation estimate of treatment effect
  N = dim(Y)[1]
  T = dim(Y)[2]
  wbar = matrix(apply(swmat,2,mean),N,T,byrow=TRUE)
  sum(Y*(swmat-wbar))/sum(swmat*(1-wbar))
}

subest = function(Y,swmat,strata){
# compute stratified permutation estimate of treatment effect
  nstrata = length(unique(strata))
  num=den=0
  for (h in 1:nstrata){
    Ystrat = Y[strata==h,]
    swstrat = swmat[strata==h,]
    N = dim(Ystrat)[1]
    T = dim(Ystrat)[2]
    wbar = matrix(apply(swstrat,2,mean),N,T,byrow=TRUE)
    num = num + sum(Ystrat*(swstrat-wbar))
    den = den + sum(swstrat*(1-wbar))
  }
  delta = num/den
  delta
}

########################
vardpc = function(Y,swmat,delta){
# corrected variance of deltahat; unbiased if var(Y_i) 
# does not depend on cluster
  N = dim(swmat)[1]
  T = dim(swmat)[2]
  wbar = apply(swmat,2,mean)
  Vpa0 = 0
  Y = Y - swmat*delta
  for (i in 1:N){
    for (j in 1:T){
      Vpa0 = Vpa0 + Y[i,j]^2*wbar[j]*(1-wbar[j])
    }
    for (j in 1:(T-1)){
      for (j1 in (j+1):T){
        Vpa0 = Vpa0 + 2*Y[i,j]*Y[i,j1]*wbar[j]*(1-wbar[j1])
      }
    }
  }
# Using Jim's way of writing sums  
#  for (i in 1:(N-1)){
#    for (i1 in (i+1):N){
#      for (j in 1:T){
#        for (j1 in 1:T) {
#          jmin = min(j,j1)
#          jmax = max(j,j1)
#          Vpa0 = Vpa0 - (2/(N-1))*Y[i,j]*Y[i1,j1]*wbar[jmin]*(1-wbar[jmax])
#        }
#      }
#    }
#  }
# Using Fan's way of writing sums (equivalent)
  for (i in 1:N){
    for (i1 in 1:N){
      if (i!=i1) {
        for (j in 1:(T-1)){
          Vpa0 = Vpa0 - (1/(N-1))*Y[i,j]*Y[i1,j]*wbar[j]*(1-wbar[j])
          for (j1 in (j+1):T){
            Vpa0 = Vpa0 - (2/(N-1))*Y[i,j]*Y[i1,j1]*wbar[j]*(1-wbar[j1])
          }
        }
          Vpa0 = Vpa0 - (1/(N-1))*Y[i,T]*Y[i1,T]*wbar[T]*(1-wbar[T])
      }
    }
  }
  norm = (N*sum(wbar*(1-wbar)))^2
  Vpa0 = Vpa0/norm
  Vpa0
}

########################
svardpc = function(Y,swmat,delta,strata){
  # corrected variance of stratified estimate of deltahat; unbiased if var(Y_i) 
  # does not depend on cluster within strata
  nstrata = length(unique(strata))
  den=Vp=0
  for (h in 1:nstrata){
    Yh = Y[strata==h,]
    N = dim(Yh)[1]
    T = dim(Yh)[2]
    wbar = apply(swmat[strata==h,],2,mean)
    Yh = Yh - swmat[strata==h,]*delta
    for (i in 1:N){
      for (j in 1:T){
        Vp = Vp + Yh[i,j]^2*wbar[j]*(1-wbar[j])
      }
      for (j in 1:(T-1)){
        for (j1 in (j+1):T){
          Vp = Vp + 2*Yh[i,j]*Yh[i,j1]*wbar[j]*(1-wbar[j1])
        }
      }
    }
    for (i in 1:N){
      for (i1 in 1:N){
        if (i!=i1) {
          for (j in 1:(T-1)){
            Vp = Vp - (1/(N-1))*Yh[i,j]*Yh[i1,j]*wbar[j]*(1-wbar[j])
            for (j1 in (j+1):T){
              Vp = Vp - (2/(N-1))*Yh[i,j]*Yh[i1,j1]*wbar[j]*(1-wbar[j1])
            }
          }
          Vp = Vp - (1/(N-1))*Yh[i,T]*Yh[i1,T]*wbar[T]*(1-wbar[T])
        }
      }
    }
    den = den + (N*sum(wbar*(1-wbar)))
  }
  svar = Vp/den^2
  svar
}


########################
vardpc3 = function(Y,swmat){
# corrected variance of deltahat; unbiased for any var(Y_i) and any delta; 
# requires at least 2 clusters with each intervention pattern
  N = dim(swmat)[1]
  T = dim(swmat)[2]
  wbar = apply(swmat,2,mean)
  #
  strata = apply(swmat,1,sum)
  svals = as.numeric(names(table(strata)))
  nnn = as.vector(table(strata))
  if (any(nnn<=1)) stop("Must have at least 2 clusters for each intervention pattern")
  #
  Vpa0 = 0
  for (k in 1:length(nnn)) {
   nstrat = nnn[k]
   istrat =  (1:N)[strata==svals[k]]
   for (ind in 1:nstrat){
    i = istrat[ind]
    for (j in 1:T){
      Vpa0 = Vpa0 + Y[i,j]^2*(swmat[i,j]-wbar[j])^2
    }
    for (j in 1:(T-1)){
      for (j1 in (j+1):T){
        Vpa0 = Vpa0 + 2*Y[i,j]*Y[i,j1]*(swmat[i,j]-wbar[j])*(swmat[i,j1]-wbar[j1])
      }
    }
  }
  for (ind in 1:(nstrat-1)){
    i = istrat[ind]
    for (ind1 in (ind+1):nstrat){
      i1 = istrat[ind1]
      for (j in 1:T){
        for (j1 in 1:T) {
          jmin = min(j,j1)
          jmax = max(j,j1)
 #         Vpa0 = Vpa0 - (2/(nstrat-1))*Y[i,j]*Y[i1,j1]*(swmat[i,jmin]-wbar[jmin])*(swmat[i1,jmax]-wbar[jmax])
          Vpa0 = Vpa0 - (2/(nstrat-1))*Y[i,j]*Y[i1,j1]*(swmat[i,j]-wbar[j])*(swmat[i1,j1]-wbar[j1])
        }
      }
    }
  }
  }
  norm = (N*sum(wbar*(1-wbar)))^2
  Vpa0 = Vpa0/norm
  Vpa0
}


########################################
confint = function(Y,swmat,level = .95){
  # Compute confidence interval based on vardpc 
  zval = abs(qnorm((1-level)/2))
  
  Zdelta = function(delta,Y,swmat,deltahat,zval){
    (deltahat - delta)/sqrt(vardpc(Y,swmat,delta)) - zval
  }
  
  del = ubest(Y,swmat)
  # used biased estimate of variance just to get a bound
  bound = 10*sqrt(vardpc(Y,swmat,del))
  upper = uniroot(Zdelta,c(del,del+bound),Y=Y,swmat=swmat,deltahat=del,zval= -zval)$root
  lower = uniroot(Zdelta,c(del-bound,del),Y=Y,swmat=swmat,deltahat=del,zval= zval)$root
  c(lower,upper)
}


gendat = function(N,T,swmat,time,delta,beta,sig,tau,eta,psi,n){
# generate cluster-period level data
  if (length(n)==1) n=rep(n,N)
  y = rep(0,N*T)
  sigma = matrix(0,N*T,N*T)
  X = model.matrix(~as.factor(time))[,-1]
  for (i in 1:N){
    i1 = (i-1)*T + 1
    i2 = (i-1)*T + T
    W = t(swmat[i,,drop=F])
    mu = delta[1] + delta[2]*W + X%*%beta
    s1 = diag(T)*sig/n[i] + matrix(tau,T,T) + eta*W%*%t(W) +psi
    y[i1:i2] = mvrnorm(1,mu,s1)
    sigma[i1:i2,i1:i2] = s1
  }    
  list(y=y,sigma=sigma)
}

gendati = function(N,T,swmat,time,delta,beta,sig,tau,eta,psi,n){
# generate individual level data  
  if (length(n)==1) n=rep(n,N)
  y = NULL
  X = model.matrix(~as.factor(time))[,-1]
  for (i in 1:N){
    W = matrix(rep(swmat[i,],rep(n[i],T)),n[i]*T,1)
    mu = delta[1] + delta[2]*W + rep(X%*%beta,rep(n[i],T))
    s1 = diag(n[i]*T)*sig/n[i] + matrix(tau,n[i]*T,n[i]*T) + eta*W%*%t(W) + psi*kronecker(diag(1,T),matrix(1,n[i],n[i]))
    y = c(y,mvrnorm(1,mu,s1))
  }    
  list(y=y)
}

gendatb = function(N,T,swmat,time,p0,delta,beta,tau,eta,psi,n){
  # generate cluster-period level data
  if (length(n)<=N) n=matrix(n,N,T)
  y = nsiz = rep(0,N*T)
  X = model.matrix(~as.factor(time))[,-1]
  a = rnorm(N,0,tau)
  b = matrix(rnorm(N*T,0,psi),N,T)
  d = rnorm(N,0,eta)
  for (i in 1:N){
    i1 = (i-1)*T + 1
    i2 = (i-1)*T + T
    W = t(swmat[i,,drop=F])
    mu = delta[1] + a[i] + b[i,] + (delta[2] + d[i])*W + X%*%beta
    y[i1:i2] = rbinom(T,n[i,],mu)
    nsiz[i1:i2] = n[i,]
  }    
  list(y=y,n=nsiz)
}

clsiz = function(N,ave,var){
  if (var==0) {
    n = rep(ave,N)
  } else {
    temp = exp(rnorm(N,log(ave),var)) + 1
    n = round(ave*temp/mean(temp))
  }
  n
}
###########################