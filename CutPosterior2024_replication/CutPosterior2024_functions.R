## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## Appendix C - functions to run CutPosterior2024_main.R ##
## Julia Kowalska, Mark van de Wiel, Stéphanie van der Pas ##

funB_sample<-function(X){
  Y=X
  X2=X[X>=0]
  X1=X[X<0]
  Y[X<0]=1/(1+exp(-2*X1))-0.5+0.4
  Y[X>=0]=(log(X2*2.+1)-0.15*X2^2)*0.6-0.20+0.4
  Y=Y+rnorm(length(X),0,0.1)
  return(Y)
}

funC_sample<-function(X){
  Y=0.48- 2.7*X+ 1.18*X^2 +1.21*X^3 + 2.54*X^4 - 3*X^5-1.9*X^6-5/(1+exp(-10*(X+1)))-10+sin(5*X-2)
  Y=Y/10+rnorm(length(X),0,0.1)
  return(Y)
}

sample_prob30<-function(X){
  P=rep(0,length(X))
  P[X>=0.]=invlogit((8.5*X[X>=0.]-1.5))/10.5+0.4-0.0007072
  P[X<0.]=(X[X<0.]+1)^4/15 +0.05 
  T=rep(0,length(X))
  for(j in 1:length(X)){
    T[j]=sample(c(1,0),1,prob=c(P[j],1-P[j]))
  }
  return(T)
}

## Bayesian models ##

# Model for the initial cutoff localization for a continuous score
#--This model fits two constant functions to the treatment data. The resulting 
# posterior samples of the cutoff c are then used to initialize LoTTA treatment model and 
# full LoTTA model.
# t - treatment data, c - cutoff, al - constant value on the left side of the cutoff, j - jump size at the cutoff
# clb - lower bound of the interval in the prior of c ~ unif(clb,cub),
# cub - upper bound of the interval in the prior of c ~ unif(clb,cub) #
cat("model
    {
  
    for ( i in 1:N ) {
      
      t[i]~dbern(paramt[i])
      paramt[i] <- ifelse(x[i]<c,al,al+j) 
    }
    
    
    al~dunif(0,1-j)
    j~dunif(jlb,1)
    c~dunif(clb,cub)
    
    
    }", file="cutoff_initial_CONT.txt") 

# Treatment model for a continuous score
#--This model fits to the treatment data two connected linear functions on each side of the cutoff. 
# It is the first stage of LoTTA.--#
# t - treatment data, c - cutoff, j - jump size at the cutoff
# a - slope, b - intercept, l/r - segment on the left/right-hand side, 1/2 - segment close/far from the cutoff 
# clb - lower bound of the interval in the prior of c ~ unif(clb,cub),
# cub - upper bound of the interval in the prior of c ~ unif(clb,cub),
# k1t - length of the window on the left-hand side, k2t - length of the window on the right-hand side
# lb - minimum window size
# ublt - minimum value of the left boundary point, ubrt - maximum value of the right boundary point
# jlb - lower bound of the jump size in the treatment probability #
cat("model
    {
  
    for ( i in 1:N ) {
      
      t[i]~dbern(paramt[i])
      paramt[i] <- ifelse(x[i]<c,ifelse(x[i]>=c-k1t,a1lt*x[i]+b1lt,a2lt*x[i]+b2lt),ifelse(x[i]<=c+k2t,a1rt*x[i]+b1rt,a2rt*x[i]+b2rt)) 
    }
    pr=0.0001
    MAX=max(x)
    MIN=min(x)
    
    
    ### Define the priors
    
    c~dunif(clb,cub)
    xc=x-c
    j~dunif(jlb,1)
    k1t~dunif(lb,c-ublt)
    k2t~dunif(lb,ubrt-c)
    a2lt~dunif(0,(1-j)/(c-k1t-MIN))
    b2lt~dunif(-a2lt*MIN,1-j-a2lt*(c-k1t))
    a1lt~dunif(0,(1-j-a2lt*(c-k1t)-b2lt)/k1t)
    b1lt=(c-k1t)*(a2lt-a1lt)+b2lt
    a1rt~dunif(0,(1-a1lt*c-b1lt-j)/k2t)
    b1rt=a1lt*c+b1lt+j-a1rt*c
    a2rt~dunif(0,(1-b1rt-(c+k2t)*a1rt)/(MAX-c-k2t))
    b2rt=(c+k2t)*(a1rt-a2rt)+b1rt
    
    
    }", file="treatment_CONT.txt") 

# LoTTA model for a continuous score and continuous outcomes
#--This model fits jointly treatment and outcome functions.--# 
# t - treatment data, y - outcome data, c - cutoff, j - jump size at the cutoff
# clb - lower bound of the interval in the prior of c ~ unif(clb,cub),
# cub - upper bound of the interval in the prior of c ~ unif(clb,cub),
# lb - minimum window size
# eff - treatment effect for compliers

# In the treatment model:
# a - slope, b - intercept, l/r - segment on the left/right-hand side, 1/2 - segment close/far from the cutoff 
# k1t - length of the window on the left-hand side, k2t - length of the window on the right-hand side
# jlb - lower bound of the jump size in the treatment probability

# In the outcome model:
# ail - coefficient of x^i on the left-hand side, air - coefficient of x^i on the right-hand side,
# kl(r) - boundary point of the window on the left(right) side of the cutoff
# ubl - minimum value of the left boundary point, ubr - maximum value of the right boundary point
# pr - precision in priors of the coefficients # 

cat("model
    {
  
    for ( i in 1:N ) {
      y[i]~dnorm(param[i],Tau[i])
      t[i]~dbern(paramt[i])
      param[i] <-ifelse(x[i]<c,a1l*(xc[i])+a0l+ilogit(100*(kl-x[i]))*((a3l)*(xc[i])^3+(a2l)*(xc[i])^2),a1r*(xc[i])+a0r+ilogit(100*(x[i]-kr))*((a3r)*(xc[i])^3+(a2r)*(xc[i])^2))
      Tau[i]<-ifelse(x[i]<c,tau1l+tau2l*ilogit(100*(kl-x[i])),tau1r+tau2r*ilogit(100*(x[i]-kr)))
      paramt[i] <- ifelse(x[i]<c,ifelse(x[i]>=c-k1t,a1lt*x[i]+b1lt,a2lt*x[i]+b2lt),ifelse(x[i]<=c+k2t,a1rt*x[i]+b1rt,a2rt*x[i]+b2rt)) 
    }
    pr=0.0001
    MAX=max(x)
    MIN=min(x)
    
    
    ### Define the priors
    c~dunif(clb,cub)
    xc=x-c
    j~dunif(jlb,1)
    k1t~dunif(lb,c-ublt)
    k2t~dunif(lb,ubrt-c)
    a2lt~dunif(0,(1-j)/(c-k1t-MIN))
    b2lt~dunif(-a2lt*MIN,1-j-a2lt*(c-k1t))
    a1lt~dunif(0,(1-j-a2lt*(c-k1t)-b2lt)/k1t)
    b1lt=(c-k1t)*(a2lt-a1lt)+b2lt
    a1rt~dunif(0,(1-a1lt*c-b1lt-j)/k2t)
    b1rt=a1lt*c+b1lt+j-a1rt*c
    a2rt~dunif(0,(1-b1rt-(c+k2t)*a1rt)/(MAX-c-k2t))
    b2rt=(c+k2t)*(a1rt-a2rt)+b1rt
    
    tau1l~dgamma(0.01,0.01)
    tau2pl~dbeta(1,1)
    tau2l=-tau2pl*tau1l
    klt~dbeta(1,1)
    kl=klt*(c-lb-ubl)+ubl
    a0l~dnorm(0,pr)
    a1l~dnorm(0,pr)
    a2l~dnorm(0,pr*(c-kl))
    a3l~dnorm(0,pr*(c-kl))
    
    tau1r~dgamma(0.01,0.01)
    tau2pr~dbeta(1,1)
    tau2r=-tau2pr*tau1r
    krt~dbeta(1,1)
    kr=krt*(ubr-c-lb)+c+lb
    a0r~dnorm(0,pr)
    a1r~dnorm(0,pr)
    a2r~dnorm(0,pr*(kr-c))
    a3r~dnorm(0,pr*(kr-c))
    
    
    
    eff=(a0r-a0l)/j
    
    }", file="LoTTA_CONT_CONT.txt") 

# LoTTA model for a sharp design and continuous outcomes
# y - outcome data, c - cutoff
# clb - lower bound of the interval in the prior of c ~ unif(clb,cub),
# cub - upper bound of the interval in the prior of c ~ unif(clb,cub),
# lb - minimum window size
# eff - treatment effect for compliers

# In the outcome model:
# ail - coefficient of x^i on the left-hand side, air - coefficient of x^i on the right-hand side,
# kl(r) - boundary point of the window on the left(right) side of the cutoff
# ubl - minimum value of the left boundary point, ubr - maximum value of the right boundary point
# pr - precision in priors of the coefficients # 

cat("model
    {
  
    for ( i in 1:N ) {
      y[i]~dnorm(param[i],Tau[i])
      param[i] <-ifelse(x[i]<c,a1l*(xc[i])+a0l+ilogit(100*(kl-x[i]))*((a3l)*(xc[i])^3+(a2l)*(xc[i])^2),a1r*(xc[i])+a0r+ilogit(100*(x[i]-kr))*((a3r)*(xc[i])^3+(a2r)*(xc[i])^2))
      Tau[i]<-ifelse(x[i]<c,tau1l+tau2l*ilogit(100*(kl-x[i])),tau1r+tau2r*ilogit(100*(x[i]-kr)))
      }
    pr=0.0001
    MAX=max(x)
    MIN=min(x)
    ### Define the priors
   
    xc=x-c
    tau1l~dgamma(0.01,0.01)
    tau2pl~dbeta(1,1)
    tau2l=-tau2pl*tau1l
    klt~dbeta(1,1)
    kl=klt*(c-lb-ubl)+ubl
    a0l~dnorm(0,pr)
    a1l~dnorm(0,pr)
    a2l~dnorm(0,pr*(c-kl))
    a3l~dnorm(0,pr*(c-kl))
    
    tau1r~dgamma(0.01,0.01)
    tau2pr~dbeta(1,1)
    tau2r=-tau2pr*tau1r
    krt~dbeta(1,1)
    kr=krt*(ubr-c-lb)+c+lb
    a0r~dnorm(0,pr)
    a1r~dnorm(0,pr)
    a2r~dnorm(0,pr*(kr-c))
    a3r~dnorm(0,pr*(kr-c))
    
    
    eff=(a0r-a0l)
    
    }", file="LoTTA_SHARP_CONT.txt")

# Function to sample an initial value of a chain

# Treatment model - continuous score
# X - score 
# T - treatment data
# C_start - posterior samples of cutoff location on a continuous scale
# obtained through "cutoff_initial_CONT.txt"
# lb - minimum window size
# ubl - minimum value of the window's left boundary point, ubr - maximum value of the window's right boundary point
# for setting initial value we recommend using ubl, ubr obtained from bounds with higher ns (ns=50)
# than for the model fitting (ns=25)
# s - seed 

Initial_treatment_CONT<-function(X,T,C_start,lb,ubr,ubl,start,prob,s,jlb=0.2){
  set.seed(s)
  pr=0.5
  MIN=min(X)
  MAX=max(X)
  c=sample(C_start,1)
  tl=mean(T[X<c])
  tr=mean(T[X>=c])
  
  
  j=max(jlb,tr-tl)
  k1t=runif(1,lb,c-ubl)
  k2t=runif(1,lb,ubr-c)
  a2lt=0
  b2lt=tl
  a1lt=0
  b1lt=(c-k1t)*(a2lt-a1lt)+b2lt
  a1rt=0
  b1rt=a1lt*c+b1lt+j-a1rt*c
  a2rt=0
  b2rt=(c+k2t)*(a1rt-a2rt)+b1rt
  return(list(c=c,j=j,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s)) 
}


# LoTTA model for a continuous score and continuous outcomes (fuzzy design)
# X - score 
# T - treatment data
# Y - outcome data
# C_start - posterior samples of cutoff location (categorized with natural numbers)
# obtained through "cutoff_initial_CONT.txt"
# ubl - minimum value of the left boundary point of the window, ubr - maximum value of the right boundary point the window
# for setting initial value we recommend using ubl, ubr obtained from bounds with higher ns (ns=50)
# than for the model fitting (ns=25)
# s - seed 
Initial_CONT_CONT<-function(X,T,Y,C_start,lb,ubr,ubl,jlb,s){
  set.seed(s)
  pr=sd(Y)
  MIN=min(X)
  MAX=max(X)
  
  c=sample(C_start,1)
  tl=mean(T[X<c])
  tr=mean(T[X>=c])
  yl=mean(Y[X<c])
  yr=mean(Y[X>=c])
  
  klt=rbeta(1,1,1)
  kl=klt*(c-lb-ubl)+ubl
  krt=rbeta(1,1,1)
  kr=krt*(ubr-c-lb)+c+lb
  kld=kl-c
  krd=kr-c
  a0l=rnorm(1,yl,pr)
  a1l=rnorm(1,0,pr)
  a2l=rnorm(1,0,pr)
  a3l=0
  a0r=rnorm(1,yr,pr)
  a1r=rnorm(1,0,pr)
  a2r=rnorm(1,0,pr)
  a3r=0
  tau1r=rchisq(1,7)
  tau2pr=rbeta(1,1,1)
  tau1l=rchisq(1,7)
  tau2pl=rbeta(1,1,1)
  j=max(tr-tl,jlb+0.0001)
  k1t=runif(1,lb+0.1*max(c-ubl-lb,0),max(c-ubl-0.1*(c-ubl-lb),lb+0.1*max(c-ubl-lb,0)+0.0001))
  k2t=runif(1,lb+0.1*max(ubr-c-lb,0),max(ubr-c-0.1*(ubr-c-lb),lb+0.1*max(ubr-c-lb,0)+0.0001))
  a2lt=0
  b2lt=tl
  a1lt=0
  b1lt=(c-k1t)*(a2lt-a1lt)+b2lt
  a1rt=0
  b1rt=a1lt*c+b1lt+j-a1rt*c
  a2rt=0
  b2rt=(c+k2t)*(a1rt-a2rt)+b1rt
  
  return(list(c=c,j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,tau1r=tau1r,tau2pr=tau2pr,tau1l=tau1l,tau2pl=tau2pl,klt=klt,krt=krt,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s)) 
}

# LoTTA model for continuous outcomes (sharp design)
# X - score 
# Y - outcome data
# c - cutoff
# lb - minimum window size
# ubl - minimum value of the left boundary point of the window, ubr - maximum value of the right boundary point the window
# for setting initial value we recommend using ubl, ubr obtained from bounds with higher ns (ns=50)
# than for the model fitting (ns=25)
# s - seed 
Initial_SHARP_CONT<-function(X,Y,c,lb,ubr,ubl,s){
  set.seed(s)
  pr=sd(Y)
  MIN=min(X)
  MAX=max(X)
  
  yl=mean(Y[X<c])
  yr=mean(Y[X>=c])
  
  klt=rbeta(1,1,1)
  kl=klt*(c-lb-ubl)+ubl
  krt=rbeta(1,1,1)
  kr=krt*(ubr-c-lb)+c+lb
  kld=kl-c
  krd=kr-c
  a0l=rnorm(1,yl,pr)
  a1l=rnorm(1,0,pr)
  a2l=rnorm(1,0,pr)
  a3l=0
  a0r=rnorm(1,yr,pr)
  a1r=rnorm(1,0,pr)
  a2r=rnorm(1,0,pr)
  a3r=0
  tau1r=rchisq(1,7)
  tau2pr=rbeta(1,1,1)
  tau1l=rchisq(1,7)
  tau2pl=rbeta(1,1,1)
  
  return(list(a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,tau1r=tau1r,tau2pr=tau2pr,tau1l=tau1l,tau2pl=tau2pl,klt=klt,krt=krt,.RNG.seed=s)) 
}

## Additional functions ##

# Function to compute minimum and maximum window size 
# X - score 
# ns - minimum number of points in the external window, default value 25
# Return list of values (ubl,ubr, lb)
# ubl - minimum value of the window's left boundary point, ubr - maximum value of the window's right boundary point
# lb - minimum window size (grid size) #
bounds<-function(X,ns=25){
  Xu=sort(unique(X))
  ql=ns/length(X)
  q=as.numeric(quantile(X,c(ql,1-ql)))
  ubr=q[2]
  ubl=q[1]
  N=length(Xu)
  diff=Xu[2:N]-Xu[1:(N-1)]
  diff1=diff[1:(N-2)]
  diff2=diff[2:(N-1)]
  Diff=ifelse(diff1>diff2,diff1,diff2)
  qd=quantile(Diff,c(0.75))
  lb=qd[[1]]
  return(list(ubl=ubl,ubr=ubr,lb=lb))
}
