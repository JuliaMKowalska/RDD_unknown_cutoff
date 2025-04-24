## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## Simulations - functions to run the Simulations_main script ##
## Julia Kowalska, Mark van de Wiel, Stéphanie van der Pas ##


## Outcome functions used in simulations ##

# Plotting

funA<-function(X){ 
  Y=X
  X2=X[X>=0]
  X1=X[X<0]
  Y[X<0]=(1.8*X1^3+2.*X1^2)+0.05
  Y[X>=0]=0.05*X2-0.1*X2^2+0.22
  return(Y)
}

funB<-function(X){
  Y=X
  X2=X[X>=0]
  X1=X[X<0]
  Y[X<0]=1/(1+exp(-2*X1))-0.5+0.4
  Y[X>=0]=(log(X2*2.+1)-0.15*X2^2)*0.6-0.20+0.4
  
  return(Y)
}

funC<-function(X){
  Y=0.48- 2.7*X+ 1.18*X^2 +1.21*X^3 + 2.54*X^4 - 3*X^5-1.9*X^6-5/(1+exp(-10*(X+1)))-10+sin(5*X-2)
  Y=Y/10
  return(Y)
}

lee<-function(X){
  X1=X[X<0]
  X2=X[X>=0]
  Y1=0.48+ 1.27*X1+ 7.18*X1^2 +20.21*X1^3 + 21.54*X1^4 + 7.33*X1^5
  Y2=0.52+ 0.84*X2 - 3.00*X2^2+7.99*X2^3 -9.01*X2^4 + 3.56*X2^5
  Y=append(Y1,Y2)
  return(Y)
}

ludwig<-function(X){
  X1=X[X<0]
  X2=X[X>=0]
  Y1=3.71+ 2.3*X1+ 3.28*X1^2 +1.45*X1^3 + 0.23*X1^4 + 0.03*X1^5
  Y2=0.26+ 18.49*X2 - 54.81*X2^2+74.3*X2^3 -45.02*X2^4 + 9.83*X2^5
  Y=append(Y1,Y2)
  return(Y) 
}

# Sampling

funA_sample<-function(X){
  Y=X
  X2=X[X>=0]
  X1=X[X<0]
  Y[X<0]=(1.8*X1^3+2.*X1^2)+0.05
  Y[X>=0]=0.05*X2-0.1*X2^2+0.22
  Y=Y+rnorm(length(X),0,0.1)
  return(Y)
}

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

lee_sample<-function(X){
  X1=X[X<0]
  X2=X[X>=0]
  Y1=0.48+ 1.27*X1+ 7.18*X1^2 +20.21*X1^3 + 21.54*X1^4 + 7.33*X1^5
  Y2=0.52+ 0.84*X2 - 3.00*X2^2+7.99*X2^3 -9.01*X2^4 + 3.56*X2^5
  Y=append(Y1,Y2)+rnorm(length(X),0,0.1295)
  return(Y)
}

ludwig_sample<-function(X){
  X1=X[X<0]
  X2=X[X>=0]
  Y1=3.71+ 2.3*X1+ 3.28*X1^2 +1.45*X1^3 + 0.23*X1^4 + 0.03*X1^5
  Y2=0.26+ 18.49*X2 - 54.81*X2^2+74.3*X2^3 -45.02*X2^4 + 9.83*X2^5
  Y=append(Y1,Y2)+rnorm(length(X),0,0.1295) 
  return(Y) 
}



## Treatment probability functions used in simulations ##
# Plotting

fun_prob55<-function(X){
  P=rep(0,length(X))
  P[X>=0.]=invlogit((8.5*X[X>=0.]-1.5))/10.5+0.65-0.0007072
  P[X<0.]=(X[X<0.]+1)^4/15+0.05 
  return(P)
}

fun_prob30<-function(X){
  P=rep(0,length(X))
  P[X>=0.]=invlogit((8.5*X[X>=0.]-1.5))/10.5+0.4-0.0007072
  P[X<0.]=(X[X<0.]+1)^4/15 +0.05 
  return(P)
}

# Sampling

sample_prob55<-function(X){
    P=rep(0,length(X))
    P[X>=0.]=invlogit((8.5*X[X>=0.]-1.5))/10.5+0.65-0.0007072
    P[X<0.]=(X[X<0.]+1)^4/15+0.05 
    T=rep(0,length(X))
    for(j in 1:length(X)){
      T[j]=sample(c(1,0),1,prob=c(P[j],1-P[j]))
    }
    return(T)
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

## Single simulation - fuzzy design
# i - seed,
# name - function name: "A", "B", "C", "lee", "ludwig",
# jump - string, jumpsize: "0.55" or "0.3", indicates treatment probablity function
simulation_FUZZY<-function(i,name,jump){                      
  functions=list("A"=funA_sample,"B"=funB_sample, "C"=funC_sample, "lee"=lee_sample, "ludwig"=ludwig_sample)
  probs=list("0.55"=sample_prob55,"0.3"=sample_prob30)
  fun=functions[[name]]
  prob=probs[[jump]]
  print(i)
  set.seed(i)
  X=sort(2*rbeta(500,2,4)-1)
  Y=fun(X)
  T=prob(X)
  (b_f1=bounds(X,25))
  b_f1t=bounds(X,25)
  ubr=b_f1$ubr
  ubl=b_f1$ubl
  ubrt=b_f1t$ubr
  ublt=b_f1t$ubl
  lb=b_f1$lb
  (b_s=bounds(X,50))
  ubrs=b_s$ubr
  ubls=b_s$ubl
  jlb=0.2
  nc=1
  dat1T=list(N=length(X),x=X,t=T,jlb=0.2,clb=-0.8,cub=0.2)
  param_full=c('c','j','kl','kr','klt','krt','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t','tau1r','tau2r','tau1l','tau2l')
  param_c=c('c')
  initc1=list(c=-0.3,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
  dat1_c<- run.jags('cutoff_initial_CONT.txt',inits = list(initc1) ,data=dat1T,monitor=param_c,burnin = 900,sample=2000,adapt = 100,n.chains = 1,method = 'simple')
  C_start=as.numeric(combine.mcmc(dat1_c$mcmc))
  init1=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,1)
  init2=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,2)
  init3=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,3)
  init4=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,4)
  dat=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,nc=nc,jlb=0.2,clb=-0.8,cub=0.2,seed=i)
  posterior<- run.jags('LoTTA_CONT_CONT.txt', data=dat,inits = list(init1,init2,init3,init4),monitor=param_full,burnin = 30000,sample=25000,adapt = 500,n.chains = 4,method = 'parallel')
}
## LLR simulations results - sharp design 
# N - number of samples,
# name - function name: "A", "B", "C", "lee", "ludwig", 
# jump - string, jumpsize: "0.55" or "0.3", indicates treatment probablity function
# st - first seed - 1 (0 corresponds to starting simulations with seed 1)
# Returns data frame with the following parameters for each simulation:
# absolute error, 95% CI length, bias,
# binary indicator if CI contains tr.eff. value, binary indicator if CI correctly identifies sign of tr.eff,
#  absolute error of compliance rate estimate #
LLR_performance_FUZZY<-function(N,name,jump,st=0){
  functions=list("A"=funA_sample,"B"=funB_sample, "C"=funC_sample, "lee"=lee_sample, "ludwig"=ludwig_sample)
  probs=list("0.55"=sample_prob55,"0.3"=sample_prob30)
  fun=functions[[name]]
  prob=probs[[jump]]
  effects=list("A"=0.17,"B"=-0.2, "C"=0, "lee"=0.04, "ludwig"=-3.45)
  compliance=list("0.55"=0.55,"0.3"=0.30)
  tr_eff=effects[[name]]/compliance[[jump]]
  rddat=data.frame(abs_err=0,ci_length=0,bias=0,cov=0,sign=0,j_abs=0)
  for(i in 1:N){
    set.seed(st+i)
    X=sort(2*rbeta(500,2,4)-1)
    Y=fun(X)
    T=prob(X)
    r=rdrobust(Y,X,0,fuzzy = T)
    meff=as.numeric(r$Estimate[2])
    merreff=abs(tr_eff-meff)
    qeff=as.numeric(r$ci[3,])
    coveff=ifelse(tr_eff<=qeff[2]&tr_eff>=qeff[1],1,0)
    if(tr_eff<0){
      signeff=ifelse(0>qeff[2],1,0)
    }
    else{
      signeff=ifelse(0<qeff[1],1,0)
    }
    qcieff_len=as.numeric(qeff[2]-qeff[1])
    j=r$tau_T[2]
    out=meff*j
    biases=tr_eff-meff
    l=list(abs_err=merreff,ci_length=qcieff_len,bias=biases,cov=coveff,sign=signeff,j_abs_err=abs(compliance[[jump]]-j))
    rddat[i,]=l
  }
  return(rddat)
}
## Results from a single posterior sample - fuzzy design 
# post_sample - posterior sample from LoTTA_CONT_CONT,
# name - function name: "A", "B", "C", "lee", "ludwig", 
# model - string, indicates Bayesian model: "LoTTA" or "3poly"  
# Returns data frame with the following parameters for each simulation:
# absolute error of median and MAP estimate of tr.eff, 95% CrI length of symmetric CrI and HDI,
# bias of median and MAP estimate,
# binary indicators if CI contains tr.eff. value for Sym.CrI and HDI, 
# binary indicator if CI correctly identifies sign of tr.eff for Sym.CrI and HDI 
# absolute error of median and MAP estimate of cutoff, 
# absolute error of median and MAP estimate of compliance rate #
performance_sample_FUZZY<-function(post_sample,name,jump){
  c=0
  functions=list("A"=funA_sample,"B"=funB_sample, "C"=funC_sample, "lee"=lee_sample, "ludwig"=ludwig_sample)
  probs=list("0.55"=sample_prob55,"0.3"=sample_prob30)
  fun=functions[[name]]
  prob=probs[[jump]]
  effects=list("A"=0.17,"B"=-0.2, "C"=0, "lee"=0.04, "ludwig"=-3.45)
  compliance=list("0.55"=0.55,"0.3"=0.30)
  j=compliance[[jump]]
  tr_eff=effects[[name]]/compliance[[jump]]
  
  Samples=combine.mcmc(post_sample)
  C=as.numeric(Samples[,1])
  J=as.numeric(Samples[,2])
  Eff=as.numeric(Samples[,7])
  
  qeff=as.numeric(quantile(Eff,c(0.025,0.975)))
  meff=median(Eff)
  mapeff=as.numeric(map_estimate(Eff))
  hdieff=as.numeric(ci(Eff,method='HDI'))[2:3]
  coveff=ifelse(tr_eff<=qeff[2]&tr_eff>=qeff[1],1,0)
  hdicoveff=ifelse(tr_eff<=hdieff[2]&tr_eff>=hdieff[1],1,0)
  
  if(tr_eff<0){
    signeff=ifelse(0>qeff[2],1,0)
    hdisigneff=ifelse(0>hdieff[2],1,0)
  }
  else{
    signeff=ifelse(0<qeff[1],1,0)
    hdisigneff=ifelse(0<hdieff[1],1,0)
  }
  
  merreff=abs(tr_eff-meff)
  maperreff=abs(tr_eff-mapeff)
  sample_err=mean(abs(tr_eff-Eff))
  qcieff_len=as.numeric(qeff[2]-qeff[1])
  hdicieff_len=as.numeric(hdieff[2]-hdieff[1])
  jes=as.numeric(map_estimate(J))
  jabs=abs(j-jes)
  ces=as.numeric(map_estimate(C))
  cabs=abs(c-ces)
  mjes=as.numeric(median(J))
  mjabs=abs(j-mjes)
  mces=as.numeric(median(C))
  mcabs=abs(c-mces)
  return(list(abs_err_med=merreff,abs_err_map=maperreff,ci_length_sym=qcieff_len,ci_length_hdi=hdicieff_len,bias_med=(tr_eff-meff),bias_map=(tr_eff-mapeff),cov_med=coveff,cov_hdi=hdicoveff,sign_med=signeff,sign_hdi=hdisigneff,c_abs_med=mcabs,c_abs_map=cabs,j_abs_med=mjabs,j_abs_map=jabs))
}
## Single simulation - sharp design LoTTA model
# i - seed,
# name - function name: "A", "B", "C", "lee", "ludwig" #
simulation_SHARP<-function(i,name){                      
  functions=list("A"=funA_sample,"B"=funB_sample, "C"=funC_sample, "lee"=lee_sample, "ludwig"=ludwig_sample)
  fun=functions[[name]]
  print(i)
  set.seed(i)
  X=sort(2*rbeta(500,2,4)-1)
  Y=fun(X)
  (b_f1=bounds(X,25))
  ubr=b_f1$ubr
  ubl=b_f1$ubl
  
  lb=b_f1$lb
  (b_s=bounds(X,50))
  ubrs=b_s$ubr
  ubls=b_s$ubl
  nc=1
  param_full_out=c('kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','tau1r','tau2r','tau1l','tau2l')
  init1=Initial_SHARP_CONT(X,Y,0,lb,ubrs,ubls,1)
  init2=Initial_SHARP_CONT(X,Y,0,lb,ubrs,ubls,2)
  init3=Initial_SHARP_CONT(X,Y,0,lb,ubrs,ubls,3)
  init4=Initial_SHARP_CONT(X,Y,0,lb,ubrs,ubls,4)
  dat2=list(N=length(X),x=X,y=Y,ubr=ubr,ubl=ubl,lb=lb,nc=nc,seed=i)
  dat12_post<- run.jags('LoTTA_SHARP_CONT.txt', data=dat2,inits = list(init1,init2,init3,init4),monitor=param_full_out,burnin = 30000,sample=25000,adapt = 500,n.chains = 4,method = 'parallel')
}

## Single simulation - sharp design cubic polynomial model
# i - seed,
# name - function name: "A", "B", "C", "lee", "ludwig" #
simulation_cubic_SHARP<-function(i,name){                      
  functions=list("A"=funA_sample,"B"=funB_sample, "C"=funC_sample, "lee"=lee_sample, "ludwig"=ludwig_sample)
  fun=functions[[name]]
  print(i)
  set.seed(i)
  X=sort(2*rbeta(500,2,4)-1)
  Y=fun(X)
  (b_f1=bounds(X,25))
  ubr=b_f1$ubr
  ubl=b_f1$ubl
  
  lb=b_f1$lb
  (b_s=bounds(X,50))
  ubrs=b_s$ubr
  ubls=b_s$ubl
  nc=1
  param_full_out=c('eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','tau1r','tau1l')
  init1=Initial_3poly_SHARP(X,Y,0,1)
  init2=Initial_3poly_SHARP(X,Y,0,2)
  init3=Initial_3poly_SHARP(X,Y,0,3)
  init4=Initial_3poly_SHARP(X,Y,0,4)
  dat2=list(N=length(X),x=X,y=Y,ubr=ubr,ubl=ubl,lb=lb,nc=nc,seed=i)
  dat12_post<- run.jags('3poly_SHARP.txt', data=dat2,inits = list(init1,init2,init3,init4),monitor=param_full_out,burnin = 30000,sample=25000,adapt = 500,n.chains = 4,method = 'parallel')
}
## Results from a single posterior sample - sharp design 
# post_sample - posterior sample from LoTTA_SHARP_CONT ot 3poly_SHARP,
# name - function name: "A", "B", "C", "lee", "ludwig", 
# model - string, indicates Bayesian model: "LoTTA" or "3poly"  
# Returns data frame with the following parameters for each simulation:
# absolute error of median and MAP estimate of tr.eff, 95% CrI length of symmetric CrI and HDI,
# bias of median and MAP estimate,
# binary indicators if CI contains tr.eff. value for Sym.CrI and HDI, 
# binary indicator if CI correctly identifies sign of tr.eff for Sym.CrI and HDI #
performance_sample_SHARP<-function(post_sample,name,model){
  effects=list("A"=0.17,"B"=-0.2, "C"=0, "lee"=0.04, "ludwig"=-3.45)
  tr_eff=effects[[name]]
  Samples=combine.mcmc(post_sample)
  if(model=='LoTTA'){
    Eff=as.numeric(Samples[,3])
  }
  if(model=='3poly'){
    Eff=as.numeric(Samples[,1])
  }
  
  qeff=as.numeric(quantile(Eff,c(0.025,0.975)))
  meff=median(Eff)
  mapeff=as.numeric(map_estimate(Eff))
  hdieff=as.numeric(ci(Eff,method='HDI'))[2:3]
  coveff=ifelse(tr_eff<=qeff[2]&tr_eff>=qeff[1],1,0)
  hdicoveff=ifelse(tr_eff<=hdieff[2]&tr_eff>=hdieff[1],1,0)
  if(tr_eff<0){
    signeff=ifelse(0>qeff[2],1,0)
    hdisigneff=ifelse(0>hdieff[2],1,0)
  }
  else{
    signeff=ifelse(0<qeff[1],1,0)
    hdisigneff=ifelse(0<hdieff[1],1,0)
  }
  
  merreff=abs(tr_eff-meff)
  maperreff=abs(tr_eff-mapeff)
  sample_err=mean(abs(tr_eff-Eff))
  qcieff_len=as.numeric(qeff[2]-qeff[1])
  hdicieff_len=as.numeric(hdieff[2]-hdieff[1])
  
  return(list(abs_err_med=merreff,abs_err_map=maperreff,ci_length_sym=qcieff_len,ci_length_hdi=hdicieff_len,bias_med=tr_eff-meff,bias_map=tr_eff-mapeff,cov_med=coveff,cov_hdi=hdicoveff,sign_med=signeff,sign_hdi=hdisigneff))
  
  
}
## LLR simulations results - sharp design 
# N - number of samples,
# name - function name: "A", "B", "C", "lee", "ludwig", 
# st - first seed - 1 (0 corresponds to starting simulations with seed 1)
# Returns data frame with the following parameters for each simulation:
# absolute error, 95% CI length, bias,
# binary indicator if CI contains tr.eff. value, binary indicator if CI correctly identifies sign of tr.eff #
LLR_performance_SHARP<-function(N,name,st=0){
  functions=list("A"=funA_sample,"B"=funB_sample, "C"=funC_sample, "lee"=lee_sample, "ludwig"=ludwig_sample)
  fun=functions[[name]]
  effects=list("A"=0.17,"B"=-0.2, "C"=0, "lee"=0.04, "ludwig"=-3.45)
  tr_eff=effects[[name]]
  rddat=data.frame(abs_err=0,ci_length=0,bias=0,cov=0,sign=0)
  
  for(i in 1:N){
    set.seed(st+i)
    X=sort(2*rbeta(500,2,4)-1)
    Y=fun(X)
    r=rdrobust(Y,X,0)
    meff=as.numeric(r$Estimate[2])
    merreff=abs(tr_eff-meff)
    qeff=as.numeric(r$ci[3,])
    coveff=ifelse(tr_eff<=qeff[2]&tr_eff>=qeff[1],1,0)
    if(tr_eff<0){
      signeff=ifelse(0>qeff[2],1,0)
    }
    else{
      signeff=ifelse(0<qeff[1],1,0)
    }
    qcieff_len=as.numeric(qeff[2]-qeff[1])
    l=list(abs_err=merreff,ci_length=qcieff_len,bias=tr_eff-meff,cov=coveff,sign=signeff)
    rddat[i,]=l
  }
  return(rddat)
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
    
    tau1l~dchisqr(7)
    tau2pl~dbeta(1,1)
    tau2l=-tau2pl*tau1l
    klt~dbeta(1,1)
    kl=klt*(c-lb-ubl)+ubl
    a0l~dnorm(0,pr)
    a1l~dnorm(0,pr)
    a2l~dnorm(0,pr*(c-kl))
    a3l~dnorm(0,pr*(c-kl))
    
    tau1r~dchisqr(7)
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
    c=0
    xc=x-c
    tau1l~dchisqr(7)
    tau2pl~dbeta(1,1)
    tau2l=-tau2pl*tau1l
    klt~dbeta(1,1)
    kl=klt*(c-lb-ubl)+ubl
    a0l~dnorm(0,pr)
    a1l~dnorm(0,pr)
    a2l~dnorm(0,pr*(c-kl))
    a3l~dnorm(0,pr*(c-kl))
    
    tau1r~dchisqr(7)
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

# Cubic polynomial model for a sharp design with continuous outcomes
# y - outcome data, c - cutoff
# clb - lower bound of the interval in the prior of c ~ unif(clb,cub),
# cub - upper bound of the interval in the prior of c ~ unif(clb,cub),
# eff - treatment effect for compliers

# In the outcome model:
# ail - coefficient of x^i on the left-hand side, air - coefficient of x^i on the right-hand side,
# pr - precision in priors of the coefficients # 

cat("model
    {
  
    for ( i in 1:N ) {
      y[i]~dnorm(param[i],Tau[i])
      param[i] <-ifelse(x[i]<c,a1l*(xc[i])+a0l+(a3l)*(xc[i])^3+(a2l)*(xc[i])^2,a1r*(xc[i])+a0r+(a3r)*(xc[i])^3+(a2r)*(xc[i])^2)
      Tau[i]<-ifelse(x[i]<c,tau1l,tau1r)
    }
    pr=0.0001
    MAX=max(x)
    MIN=min(x)
    
    
    ### Define the priors
    c=0
    xc=x-c
    
    a0l~dnorm(0,pr)
    a1l~dnorm(0,pr)
    a2l~dnorm(0,pr)
    a3l~dnorm(0,pr)
    tau1l~dchisqr(7)
    
    a0r~dnorm(0,pr)
    a1r~dnorm(0,pr)
    a2r~dnorm(0,pr)
    a3r~dnorm(0,pr)
    tau1r~dchisqr(7)
    
    eff=(a0r-a0l)
    
    }", file="3poly_SHARP.txt") 
# Function to sample an initial value of a chain

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
# Cubic model for continuous outcomes (sharp design)
# X - score 
# Y - outcome data
# c - cutoff
# s - seed 
Initial_3poly_SHARP<-function(X,Y,c,s){
  set.seed(s)
  pr=sd(Y)
  MIN=min(X)
  MAX=max(X)
  
  yl=mean(Y[X<c])
  yr=mean(Y[X>=c])
  
  a0l=rnorm(1,yl,pr)
  a1l=rnorm(1,0,pr)
  a2l=rnorm(1,0,pr)
  a3l=0
  a0r=rnorm(1,yr,pr)
  a1r=rnorm(1,0,pr)
  a2r=rnorm(1,0,pr)
  a3r=0
  tau1r=rchisq(1,7)
  tau1l=rchisq(1,7)
  
  return(list(a0l=a0l,a1l=a1l,a2l=a2l,a3l=a3l,a0r=a0r,a1r=a1r,a2r=a2r,a3r=a3r,tau1r=tau1r,tau1l=tau1l,.RNG.seed=s)) 
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

# Logit and inverse logit functions
logit<-function(x){
  return(log(x/(1-x)))
}
invlogit<-function(x){
  return(1/(1+exp(-x)))
}

