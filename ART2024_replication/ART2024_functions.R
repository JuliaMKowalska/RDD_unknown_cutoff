## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## ART Application - functions to run ART2024_main.R and ART2024_cont.R ##
## Julia Kowalska, Mark van de Wiel, Stéphanie van der Pas ##

## Bayesian models ##

# Model for the initial cutoff localization for a discrete score
#--This model fits two constant functions to the treatment data. The resulting 
# posterior samples of the cutoff c are then used to initialize LoTTA treatment model and 
# full LoTTA model.
# t - treatment data, c - cutoff, al - constant value on the left side of the cutoff, j - jump size at the cutoff
# prob - vector of length n of prior weights on the discrete location of c,
# ct - categorical variable, with discrete distribution accoridng to prob on 1,...,n
# start, nc - shifting and rescaling coefficient to translate prior on ct to the domain of c--#
cat("model
    {
  
    for ( i in 1:N ) {
      
      t[i]~dbern(paramt[i])
      paramt[i] <- ifelse(x[i]<c,al,al+j) 
    }
    
    
    al~dunif(0,1-j)
    j~dunif(jlb,1)
    c=(ct+start)/nc
    ct~dcat(prob)
    
    
    }", file="cutoff_initial_DIS.txt") 

# Treatment model for a discrete score
#--This model fits to the treatment data two connected linear functions on each side of the cutoff. 
# It is the first stage of LoTTA.--#
# t - treatment data, c - cutoff, j - jump size at the cutoff
# a - slope, b - intercept, l/r - segment on the left/right-hand side, 1/2 - segment close/far from the cutoff 
# k1t - length of the window on the left-hand side, k2t - length of the window on the right-hand side
# lb - grid size
# ublt - minimum value of the left boundary point, ubrt - maximum value of the right boundary point
# prob - vector of length n of prior weights on the discrete location of c,
# ct - categorical variable, with discrete distribution according to prob on 1,...,n
# start, nc - shifting and rescaling coefficient to translate prior on ct to the domain of c 
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
    ct~dcat(prob)
    c=(ct+start)/nc
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
    
    
    }", file="treatment_DIS.txt") 
# LoTTA model for a discrete score and binary outcomes
#--This model fits jointly treatment and outcome functions.--# 
# t - treatment data, y - outcome data, c - cutoff, j - jump size at the cutoff
# prob - vector of length n of prior weights on the discrete location of c,
# ct - categorical variable, with discrete distribution accoridng to prob on 1,...,n
# start, nc - shifting and rescaling coefficient to translate prior on ct to the domain of c--#
# lb - grid size
# eff - treatment effect for compliers

# In the treatment model:
# a - slope, b - intercept, l/r - segment on the left/right-hand side, 1/2 - segment close/far from the cutoff 
# k1t - length of the window on the left-hand side, k2t - length of the window on the right-hand side
# jlb - lower bound of the jump size in the treatment probability

# In the outcome model:
# ail - coefficient of x^i on the left-hand side, air - coefficient of x^i on the right-hand side,
# b0l(r), b1l(r) - coefficients of the linear part in the logit function on the left(right)-hand side
# calculated so that a0l(r)+a1l(r)*(x-c) is the first order approximation of the function in the tail
# on the left(right) side 
# kl(r) - boundary point of the window on the left(right) side of the cutoff
# ubl - minimum value of the left boundary point, ubr - maximum value of the right boundary point
# pr - precision in priors of the coefficients # 

cat("model
    {
  
    for ( i in 1:N ) {
      y[i]~dbern(param[i])
      t[i]~dbern(paramt[i])
      param[i] <-ifelse(x[i]<c,a1l*(xc[i])+a0l+ilogit(100*(kl-x[i]))*(ilogit((a3l)*(xc[i])^3+(a2l)*(xc[i])^2+b1l*(xc[i])+b0l)-a0l-a1l*xc[i]),a1r*(xc[i])+a0r+ilogit(100*(x[i]-kr))*(ilogit((a3r)*(xc[i])^3+(a2r)*(xc[i])^2+b1r*(xc[i])+b0r)-a0r-a1r*xc[i]))
      paramt[i] <- ifelse(x[i]<c,ifelse(x[i]>=c-k1t,a1lt*x[i]+b1lt,a2lt*x[i]+b2lt),ifelse(x[i]<=c+k2t,a1rt*x[i]+b1rt,a2rt*x[i]+b2rt)) 

    }
    pr=0.00001
    MAX=max(x)
    MIN=min(x)
    ct~dcat(prob)
    c=(ct+start)/nc
    xc=x-c
    j~dunif(max(jlb,abs(a0r-a0l)),1)
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
    
    kl~dunif(ubl,c-lb)
    kld=kl-c
    a0l~dunif(0,1)
    b0l=logit(a0l)
    a1l~dunif((0.99-a0l)/(kld),(0.01-a0l)/(kld))
    b1l=a1l*pow(ilogit(-b0l)*(1-ilogit(-b0l)),-1)
    a2l~dnorm(0,-pr*kld)
    a3l~dnorm(0,-pr*kld)
    
    kr~dunif(c+lb,ubr)
    krd=kr-c
    a0r~dunif(0,1)
    b0r=logit(a0r)
    a1r~dunif((0.01-a0r)/(krd),(0.99-a0r)/(krd))
    b1r=a1r*pow(ilogit(-b0r)*(1-ilogit(-b0r)),-1)
    a2r~dnorm(0,pr*krd)
    a3r~dnorm(0,pr*krd)
    eff=(a0r-a0l)/j
    
    
    }", file="LoTTA_DIS_BIN.txt")

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
# Function to sample an initial value of a chain
# Treatment model disecrete score
# X - score 
# T - treatment data
# Ct_start - posterior samples of cutoff location (categorized with natural numbers)
# obtained through "cutoff_initial_dis.txt"
# ubl - minimum value of the window's left boundary point, ubr - maximum value of the window's right boundary point
# for setting initial value we recommend using ubl, ubr obtained from bounds with higher ns (ns=50)
# than for the model fitting (ns=25)
# s - seed 

Initial_treatment_DIS<-function(X,T,Ct_start,lb,ubr,ubl,start,prob,s,nc,jlb=0.2){
  set.seed(s)
  pr=0.5
  MIN=min(X)
  MAX=max(X)
  ct=sample(Ct_start,1)
  c=(ct+start)/nc
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
  return(list(ct=ct,j=j,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s)) 
}

# LoTTA model for a discrete score and binary outcomes
# X - score 
# T - treatment data
# Y - outcome data
# Ct_start - posterior samples of cutoff location (categorized with natural numbers)
# obtained through "cutoff_initial_DIS.txt"
# ubl - minimum value of the left boundary point of the window, ubr - maximum value of the right boundary point the window
# for setting initial value we recommend using ubl, ubr obtained from bounds with higher ns (ns=50)
# than for the model fitting (ns=25)
# s - seed 
Initial_DIS_BIN<-function(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,s,nc,jlb=0.2){
  set.seed(s)
  pr=0.5
  MIN=min(X)
  MAX=max(X)
  ct=sample(Ct_start,1)
  c=(ct+start)/nc
  tl=mean(T[X<c])
  tr=mean(T[X>=c])
  yl=mean(Y[X<c])
  yr=mean(Y[X>=c])
  
  kl=runif(1,ubl,c-lb)
  kr=runif(1,c+lb,ubr)
  krd=kr-c
  kld=kl-c
  a0l=runif(1,0.9*yl,min(1,1.1*yl))
  a1l=runif(1,(yl-a0l)/kld,(0.9*yl-a0l)/kld)
  b0l=logit(a0l)
  b1l=a1l*1/(invlogit(b0l)*(1-invlogit(b0l)))
  a2l=rnorm(1,0,pr)
  a3l=0
  a0r=runif(1,0.9*yr,min(1,1.1*yr))
  b0r=logit(a0r)
  a1r=runif(1,(0.9*yr-a0r)/krd,(yr-a0r)/(kr-c))
  b1r=a1r*1/(invlogit(b0r)*(1-invlogit(b0r)))
  a2r=rnorm(1,0,pr)
  a3r=0
  
  
  j=ifelse(tr-tl>=abs(a0r-a0l),max(jlb,tr-tl),max(jlb,abs(a0r-a0l)))
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
  return(list(ct=ct,j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,kl=kl,kr=kr,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s)) 
}


## Functions to plot the data ##
# Function to bin the data
#-- Outcome Y is binned with respect to the binning in the domain X. --#
# Y - data on y-axis
# X - data on x-axis
# c - cutoff, if not NULL the data is binned seperately on the two sides of the cutoff
# binsize - length of intervals in which data is binned
# The function returns a list ('Y_binned','X_binned')
# Y_binned: binned Y, X_binned - binned X
Bin_data<-function(Y,X,c=NULL,binsize=0.2){
  if(is.null(c)==1){
    L=seq(min(X),max(X),binsize)
  }
  else{
    L=seq(min(X),c,binsize)
    L=append(L,seq(c,max(X),binsize))
  }
  
  L=append(L,c(max(X)+0.0000001))
  L=unique(L)
  groups <- cut(X, L,right=FALSE,labels = L[-length(L)])
  Group<-as.numeric(groups)
  
  Binned=lapply(split(Y, groups), mean)
  Binned=as.numeric(Binned)
  return(list('Y_binned'=Binned,'X_binned'= L[-length(L)]))
}
# Logit and inverse logit functions
logit<-function(x){
  return(log(x/(1-x)))
}
invlogit<-function(x){
  return(1/(1+exp(-x)))
}
# Function that returns values of a posterior outcome function in binary LoTTA
# coef_s - vector, one posterior sample
# x - domain of the posterior function
# nc - rescaling parameter, nc=1 if the score X was not rescaled
binary_outcome_function_sample<-function(coef_s,x,nc=1){
  x=x/nc
  c=coef_s['c']
  a0l=coef_s['a0l']
  a1l=coef_s['a1l']
  a2l=coef_s['a2l']
  a3l=coef_s['a3l']
  
  a0r=coef_s['a0r']
  a1r=coef_s['a1r']
  a2r=coef_s['a2r']
  a3r=coef_s['a3r']
  
  
  b0l=logit(a0l)
  b1l=a1l*(invlogit(-b0l)*(1-invlogit(-b0l)))^(-1)
  
  b0r=logit(a0r)
  b1r=a1r*(invlogit(-b0r)*(1-invlogit(-b0r)))^(-1)
  
  kl=coef_s['kl']
  kr=coef_s['kr']
  xc=x-c
  param <-ifelse(x<c,a1l*(xc)+a0l+invlogit(100*(kl-x))*(invlogit((a3l)*(xc)^3+(a2l)*(xc)^2+b1l*(xc)+b0l)-a0l-a1l*xc),a1r*(xc)+a0r+invlogit(100*(x-kr))*(invlogit((a3r)*(xc)^3+(a2r)*(xc)^2+b1r*(xc)+b0r)-a0r-a1r*xc))
  return(param)
}
# Function that returns values of a treatment probability function in LoTTA 
# coef_s - vector, one posterior sample
# x - domain of the posterior function
# nc - rescaling parameter, nc=1 if the score X was not rescaled
treatment_function_sample<-function(coef_s,x,nc){
  x=x/nc
  c=coef_s['c']
  
  a1lt=coef_s['a1lt']
  a2lt=coef_s['a2lt']
  
  a1rt=coef_s['a1rt']
  a2rt=coef_s['a2rt']
  
  b1lt=coef_s['b1lt']
  b2lt=coef_s['b2lt']
  
  b1rt=coef_s['b1rt']
  b2rt=coef_s['b2rt']
  
  k1t=coef_s['k1t']
  k2t=coef_s['k2t']
  
  paramt<- ifelse(x<c,ifelse(x>=c-k1t,a1lt*x+b1lt,a2lt*x+b2lt),ifelse(x<=c+k2t,a1rt*x+b1rt,a2rt*x+b2rt)) 
  
  return(paramt)
}

## Functions for ART2024_cont.R ##
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

# LoTTA model for a continuous score and binary outcomes
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
# b0l(r), b1l(r) - coefficients of the linear part in the logit function on the left(right)-hand side
# calculated so that a0l(r)+a1l(r)*(x-c) is the first order approximation of the function in the tail
# on the left(right) side 
# kl(r) - boundary point of the window on the left(right) side of the cutoff
# ubl - minimum value of the left boundary point, ubr - maximum value of the right boundary point
# pr - precision in priors of the coefficients # 

cat("model
    {
  
    for ( i in 1:N ) {
      y[i]~dbern(param[i])
      t[i]~dbern(paramt[i])
      param[i] <-ifelse(x[i]<c,a1l*(xc[i])+a0l+ilogit(100*(kl-x[i]))*(ilogit((a3l)*(xc[i])^3+(a2l)*(xc[i])^2+b1l*(xc[i])+b0l)-a0l-a1l*xc[i]),a1r*(xc[i])+a0r+ilogit(100*(x[i]-kr))*(ilogit((a3r)*(xc[i])^3+(a2r)*(xc[i])^2+b1r*(xc[i])+b0r)-a0r-a1r*xc[i]))
      paramt[i] <- ifelse(x[i]<c,ifelse(x[i]>=c-k1t,a1lt*x[i]+b1lt,a2lt*x[i]+b2lt),ifelse(x[i]<=c+k2t,a1rt*x[i]+b1rt,a2rt*x[i]+b2rt)) 

    }
    pr=0.00001
    MAX=max(x)
    MIN=min(x)
    c~dunif(clb,cub)
    xc=x-c
    j~dunif(max(jlb,abs(a0r-a0l)),1)
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
    
    kl~dunif(ubl,c-lb)
    kld=kl-c
    a0l~dunif(0,1)
    b0l=logit(a0l)
    a1l~dunif((0.99-a0l)/(kld),(0.01-a0l)/(kld))
    b1l=a1l*pow(ilogit(-b0l)*(1-ilogit(-b0l)),-1)
    a2l~dnorm(0,-pr*kld)
    a3l~dnorm(0,-pr*kld)
    
    kr~dunif(c+lb,ubr)
    krd=kr-c
    a0r~dunif(0,1)
    b0r=logit(a0r)
    a1r~dunif((0.01-a0r)/(krd),(0.99-a0r)/(krd))
    b1r=a1r*pow(ilogit(-b0r)*(1-ilogit(-b0r)),-1)
    a2r~dnorm(0,pr*krd)
    a3r~dnorm(0,pr*krd)
    eff=(a0r-a0l)/j
    
    
    }", file="LoTTA_CONT_BIN.txt")

# Functions to sample an initial value of a chain
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

Initial_treatment_CONT<-function(X,T,C_start,lb,ubr,ubl,s,jlb=0.2){
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

# LoTTA model 
# X - score 
# T - treatment data
# Y - outcome data
# C_start - posterior samples of cutoff location on a continuous scale
# obtained through "cutoff_initial_dis.txt"
# ubl - minimum value of the left boundary point of the window, ubr - maximum value of the right boundary point the window
# for setting initial value we recommend using ubl, ubr obtained from bounds with higher ns (ns=50)
# than for the model fitting (ns=25)
# s - seed 
Initial_CONT_BIN<-function(X,T,Y,C_start,lb,ubr,ubl,s,jlb=0.2){
  set.seed(s)
  pr=0.5
  MIN=min(X)
  MAX=max(X)
  c=sample(C_start,1)
  tl=mean(T[X<c])
  tr=mean(T[X>=c])
  yl=mean(Y[X<c])
  yr=mean(Y[X>=c])
  
  kl=runif(1,ubl,c-lb)
  kr=runif(1,c+lb,ubr)
  krd=kr-c
  kld=kl-c
  a0l=runif(1,0.9*yl,min(1,1.1*yl))
  a1l=runif(1,(yl-a0l)/kld,(0.9*yl-a0l)/kld)
  b0l=logit(a0l)
  b1l=a1l*1/(invlogit(b0l)*(1-invlogit(b0l)))
  a2l=rnorm(1,0,pr)
  a3l=0
  a0r=runif(1,0.9*yr,min(1,1.1*yr))
  b0r=logit(a0r)
  a1r=runif(1,(0.9*yr-a0r)/krd,(yr-a0r)/(kr-c))
  b1r=a1r*1/(invlogit(b0r)*(1-invlogit(b0r)))
  a2r=rnorm(1,0,pr)
  a3r=0
  
  
  j=ifelse(tr-tl>=abs(a0r-a0l),max(jlb,tr-tl),max(jlb,abs(a0r-a0l)))
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
  return(list(c=c,j=j,a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,kl=kl,kr=kr,k1t=k1t,k2t=k2t,a1lt=a1lt,a2lt=a2lt,b2lt=b2lt,a1rt=a1rt,a2rt=a2rt,.RNG.seed=s)) 
}

# Appendix E2

# Function to calculate Y(Z=1,X=x)-Y(Z=0,X=x) from a single posterior sample, binary LoTTA 
# coef_s - posterior sample from the LoTTA model
# x - list of points 
# nc - normalization constant; default value 1
# Returns list of numerical values
binary_outcome_extrapolation<-function(coef_s,x,nc=1){
  x=x/nc
  c=coef_s['c']
  a0l=coef_s['a0l']
  a1l=coef_s['a1l']
  a2l=coef_s['a2l']
  a3l=coef_s['a3l']
  
  a0r=coef_s['a0r']
  a1r=coef_s['a1r']
  a2r=coef_s['a2r']
  a3r=coef_s['a3r']
  
  b0l=logit(a0l)
  b1l=a1l*(invlogit(-b0l)*(1-invlogit(-b0l)))^(-1)
  
  b0r=logit(a0r)
  b1r=a1r*(invlogit(-b0r)*(1-invlogit(-b0r)))^(-1)
  
  kl=coef_s['kl']
  kr=coef_s['kr']
  xc=x-c
  f_l=a1l*(xc)+a0l+invlogit(100*(kl-x))*(invlogit((a3l)*(xc)^3+(a2l)*(xc)^2+b1l*(xc)+b0l)-a0l-a1l*xc)
  f_r= a1r*(xc)+a0r+invlogit(100*(x-kr))*(invlogit((a3r)*(xc)^3+(a2r)*(xc)^2+b1r*(xc)+b0r)-a0r-a1r*xc)
  return(f_r-f_l)
}

# Function to calculate T(Z=1,X=x)-T(Z=0,X=x) from a single posterior sample 
# coef_s - posterior sample from the LoTTA model
# x - list of points 
# nc - normalization constant; default value 1
# Returns list of numerical values
treatment_function_extrapolation<-function(coef_s,x,nc=1){
  x=x/nc
  c=coef_s['c']
  
  a1lt=coef_s['a1lt']
  a2lt=coef_s['a2lt']
  
  a1rt=coef_s['a1rt']
  a2rt=coef_s['a2rt']
  
  b1lt=coef_s['b1lt']
  b2lt=coef_s['b2lt']
  
  b1rt=coef_s['b1rt']
  b2rt=coef_s['b2rt']
  
  k1t=coef_s['k1t']
  k2t=coef_s['k2t']
  
  t_l=ifelse(x>=c-k1t,a1lt*x+b1lt,a2lt*x+b2lt)
  t_r= ifelse(x<=c+k2t,a1rt*x+b1rt,a2rt*x+b2rt)
  return(t_r-t_l)
}

# Function to calculate the extrapolated treatment effect in the ART 
# Samples - combined posterior samples from the LoTTA model fitted to ART dataset
# h - window size
# data - ART dataset
# Returns list of the posterior extrapolated treatment effect 
treatment_effect_window<-function(Samples,h,data){
  tr_eff=c()
  X=data$cd4/1000
  h=h/1000
  for(i in 1:nrow(Samples)){
   
    c=Samples[i,'c']
    x=X[(X-c)<h&(c-X)<=h]
    Tef=treatment_function_extrapolation(Samples[i,],x)
    Yef=binary_outcome_extrapolation(Samples[i,],x)
    tr_eff[i]=mean(Yef)/mean(Tef)
  }
  return(tr_eff)
}

