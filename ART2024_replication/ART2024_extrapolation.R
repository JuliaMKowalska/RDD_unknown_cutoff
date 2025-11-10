## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## ART Application - Appendix E2 ##
## Julia Kowalska, Mark van de Wiel, Stéphanie van der Pas ##

## To run the following code it is necessary to first install JAGS. The instructions
# can be found under the link https://mcmc-jags.sourceforge.io ##

#install.packages('foreign')
#install.packages('ggplot2')
#install.packages('ggpubr')
#install.packages('runjags')
#install.packages('rjags')
#install.packages('bayestestR')
#install.packages('dplyr')
#install.packages('rdrobust')

library(foreign)
library(ggplot2)
library(ggpubr)
library(runjags)
library(rjags)
library(bayestestR)
library(dplyr)
library(rdrobust)

# load the required functions 

source('ART2024_functions.R')

## The data was made available by Matias D. Cattaneo and can be downloaded under the link:
# https://github.com/rdpackages-replication/CKT_2023_SIM/blob/master/CKT_2023_SIM--ART.dta ##

#download.file("https://github.com/rdpackages-replication/CKT_2023_SIM/blob/master/CKT_2023_SIM--ART.dta",destfile = "CKT_2023_SIM--ART.dta")
data <- read.dta("CKT_2023_SIM--ART.dta")
dd = data[complete.cases(data[,c("visit_test_6_18", "art_6m", "cd4")]),]
X=dd$cd4
Y=dd$visit_test_6_18
T=dd$art_6m
T=ifelse(T==1,0,1)

## Main code ##
# Step 1: Running LoTTA model 
nc=1000 # normalizing constant
ddtr=subset(dd,dd$cd4>=50 & dd$cd4<=950)
X=ddtr$cd4/nc
Y=ddtr$visit_test_6_18
T=ddtr$art_6m
T=ifelse(T==1,0,1)

# Window bounds
{b_hiv=bounds(X,25)
  ubr=b_hiv$ubr
  ubl=b_hiv$ubl
  lb=b_hiv$lb
  b_hivt=bounds(X,25)
  ubrt=b_hivt$ubr
  ublt=b_hivt$ubl}

# Prior on c - uniform between 300 and 399
start=299
prob=rep(1,100) 
# Lower bound on the jump size in the treatment probability function
jlb=0.2

# Fitting two constant functions to initialize values of c (approx 50 seconds)
datHIV_T=list(N=length(X),x=X,t=T,jlb=0.2,prob=prob,start=start,nc=nc,lb=lb)
initcART=list(ct=50,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
param_c=c('ct')
datHIV_c<- run.jags('cutoff_initial_DIS.txt',inits = list(initcART) ,data=datHIV_T,monitor=param_c,burnin = 1000,sample=2000,adapt = 100,n.chains = 1,method = 'simple')
Ct_start=as.numeric(combine.mcmc(datHIV_c))

# Fitting treatment model (approx 35 minutes)
init1=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,100,nc)
init2=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,200,nc)
init3=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,300,nc)
init4=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,400,nc)
param_cj=c('c','j')
datHIV_T=list(N=length(X),x=X,t=T,jlb=jlb,prob=prob,start=start,nc=nc,lb=lb,ubrt=ubrt,ublt=ublt)
system.time(HIV_FULL_treatment<- run.jags('treatment_DIS.txt', data=datHIV_T,monitor=param_cj,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))

# Convergence check
plot(HIV_FULL_treatment,c('trace'))

# Fitting full LoTTA model (approx 13.5 hours)
init1=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,100,nc)
init2=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,200,nc)
init3=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,300,nc)
init4=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,400,nc)
datHivfull=list(N=length(X),x=X,t=T,y=Y,ubr=b_hiv$ubr,ubl=b_hiv$ubl,ubrt=ubrt,ublt=ublt,lb=b_hiv$lb,prob=prob,start=start,nc=nc,jlb=jlb)
param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
system.time(HIV_FULL<- run.jags('LoTTA_DIS_BIN', data=datHivfull,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))
# Convergence check
plot(HIV_FULL,c('trace'),c('eff','c','j'))
gelman.diag(HIV_FULL)

# Step 2: Covariate check in the chosen windows h = 4,10,20 around MAP cutoff = 355

# The model for the covariate check, b[i] is the probability that the averages of the 
# ith covariate are the same on the both sides of the cutoff,
# for the details see Li, F., Mattei, A., and Mealli, F. (2015).

cat("model
    {
  
    for ( i in 1:N ) {
      cov[i]~dbern(paramcov[i])
      
      paramcov[i] <- 1 - pnorm(0, mu1[group[i]]+mu2[group[i]]*z[i], 1)}
    for (g in 1:G){
      
      b[g]~dbern(p)
      mu1[g] ~ dnorm(0, 0.001)
      mu2_tilde[g] ~ dnorm(0, 0.001) 
      mu2[g]=(1-b[g])*mu2_tilde[g] 
    }
    p~dbeta(1,1)
    pr=0.001
    mu_pop1~dnorm(0,pr)
    mu_pop2~dnorm(0,pr)
    tau1~dgamma(0.01,0.01)
    tau2~dgamma(0.01,0.01)
    
    for(g in 1:G){
    eff1[g] <- 1 - pnorm(0, mu1[g], 1)
    eff2[g] <- 1 - pnorm(0, mu1[g]+mu2[g], 1)
    }
    
    
    }", file="balance_binarycov.txt") 

cov=c("age1", "age2", "age3", "age4",
      "age5", "age6", "age7", "age8",
      "qtr1", "qtr2", "qtr3", "qtr4",
      "qtr5", "qtr6", "female", "clinic_a",
      "clinic_b", "clinic_c") # list of covariates


h=4
#h=10
#h=20
ddh=subset.data.frame(dd, cd4<355+h & cd4>=355-h)
ddh[,'z']=ifelse(ddh$cd4>=355,1,0) # z=1 if Xi>=c and z=0 if Xi<c;
Cov=c()
Z=c()
group=c()
for(i in 1:length(cov)){
  cov1=ddh[,c(cov[i],'z')]
  Cov1=as.numeric(as.matrix(cov1[,cov[i]]))
  Z1=as.numeric(as.matrix(cov1[,'z']))
  Cov=append(Cov,Cov1)
  Z=append(Z,Z1)
  group=append(group,rep(i,length(Z1)))
}

dat1Cov=list(N=length(Cov),cov=Cov,G=length(cov),group=group,z=Z)
param_b=c('b')
posterior<- run.jags('balance_binarycov.txt',data=dat1Cov,monitor=param_b,burnin = 10000,sample=25000,adapt = 1000,n.chains = 4,method = 'parallel')
posterior
  

# Step 3: Check the quantiles of the LoTTA windows

Samples=combine.mcmc(HIV_FULL)
kr=(Samples[,'kr']-Samples[,'c'])*1000
kl=(Samples[,'c']-Samples[,'kl'])*1000
kr=(Samples[,'kr'])*1000
kl=(Samples[,'c']-Samples[,'kl'])*1000
krt=(Samples[,'k2t'])*1000
klt=(Samples[,'k1t'])*1000
quantile(kl,c(0.05)) # above 20
quantile(kr,c(0.05)) # above 20
quantile(klt,c(0.05)) # above 4
quantile(krt,c(0.05)) # above 20
?quantile
# quantile position of h = 10 and h = 20
mean(ifelse(klt<10,1,0))
mean(ifelse(klt<20,1,0))

# Step 4: extrapolate the treatment effect to h = 4,10,20

for(h in c(4,10,20)){
  eff=treatment_effect_window(Samples,h,dd)
  mapeff=as.numeric(map_estimate(eff))
  hdieff=as.numeric(ci(eff,method='HDI'))[2:3]
  print(paste("h =",as.character(h)))
  print(paste("MAP treatment effect =",as.character(round(mapeff,3)), ' Credible interval: (',round(hdieff[1],3),',',round(hdieff[2],3),")"))
}
