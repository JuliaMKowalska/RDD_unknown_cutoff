## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## ART Application - main code ##
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

source('CutPosterior2024_functions.R')

## Negative treatment effect
# Sample the data

i=5000
set.seed(i)
X=sort(2*rbeta(500,2,2)-1)
Y=funB_sample(X)
T=sample_prob30(X)
{b_f1=bounds(X,25)
  b_f1t=bounds(X,25)
  ubr=b_f1$ubr
  ubl=b_f1$ubl
  ubrt=b_f1t$ubr
  ublt=b_f1t$ubl
  lb=b_f1$lb
  b_s=bounds(X,50)
  ubrs=b_s$ubr
  ubls=b_s$ubl}
jlb=0.2
nc=1
# Set prior on c - uniform between -0.8 and 0.2
clb=-0.8
cub=0.2

dat1T=list(N=length(X),x=X,t=T,jlb=jlb,clb=-0.8,cub=0.2)
param_full=c('c','j','kl','kr','klt','krt','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t','tau1r','tau2r','tau1l','tau2l')
param_c=c('c')
initc1=list(c=-0.3,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
dat1_c<- run.jags('cutoff_initial_CONT.txt',inits = list(initc1) ,data=dat1T,monitor=param_c,burnin = 900,sample=2000,adapt = 100,n.chains = 1,method = 'simple')
C_start=as.numeric(combine.mcmc(dat1_c$mcmc))

#Fit treatment model
init1=Initial_treatment_CONT(X,T,C_start,lb,ubr,ubl,start,prob,100)
init2=Initial_treatment_CONT(X,T,C_start,lb,ubr,ubl,start,prob,200)
init3=Initial_treatment_CONT(X,T,C_start,lb,ubr,ubl,start,prob,300)
init4=Initial_treatment_CONT(X,T,C_start,lb,ubr,ubl,start,prob,400)

dat2=list(N=length(X),x=X,t=T,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,nc=nc,jlb=jlb,clb=clb,cub=cub)
dattreatment20<-run.jags('treatment_CONT.txt', data=dat2,inits = list(init1,init2,init3,init4),monitor=param_full,burnin = 30000,sample=250,adapt = 500,n.chains = 4,method = 'parallel')

# Fit joint model
init1=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,100)  
init2=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,200)
init3=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,300)
init4=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,400)
dat1=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,nc=nc,jlb=jlb,clb=clb,cub=cub)
param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
system.time(datjoint20<- run.jags('LoTTA_CONT_CONT', data=dat1,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 30000,sample=250,adapt = 500, method='parallel',n.chains = 4))

# Sample from cut posterior
Samples=combine.mcmc(dattreatment20)
C=as.numeric(Samples[,1])
J=as.numeric(Samples[,2])
Eff_joint=as.numeric(Samples[,7])
hist(Eff_joint)
Eff_cutposterior20=c()

for(i in 1:1000){
  print(i)
  c=C[i]
  init1=Initial_SHARP_CONT(X,Y,c,lb,ubrs,ubls,1)
  init2=Initial_SHARP_CONT(X,Y,c,lb,ubrs,ubls,2)
  init3=Initial_SHARP_CONT(X,Y,c,lb,ubrs,ubls,3)
  init4=Initial_SHARP_CONT(X,Y,c,lb,ubrs,ubls,4)
  dat2=list(N=length(X),x=X,y=Y,c=c,ubr=ubr,ubl=ubl,lb=lb)
  param=c('eff')
  datout<-run.jags('LoTTA_SHARP_CONT.txt', data=dat2,inits = list(init1,init2,init3,init4),monitor=param,burnin = 10000,sample=1,adapt = 500,n.chains = 1,method = 'parallel')
  Sampleseff=combine.mcmc(datout)
  eff=as.numeric(Sampleseff[1])/J[i]
  Eff_cutposterior20=append(Eff_cutposterior20,eff)
}
# Plot the results
Samples_joint=combine.mcmc(datjoint20)
C_joint20=as.numeric(Samples_joint[,1])
c_j=data.frame(C=C_joint20)
c_j$Model<-'Joint'

c_c=data.frame(C=C)
c_c$Model<-'Cut posterior'

count.df <- rbind(c_j,c_c)
plot <- ggplot(count.df, aes(C,fill=Model))
(pl1=plot + geom_histogram(aes(y=..density..),position = 'identity') +xlim(-0.08,max(append(C_joint20,C)))+ 
  scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5)))+geom_vline(xintercept = 0,linetype='dotted',color='black') +labs(title="Posterior probability of the cutoff location ",y="density",
                                                                                                                                         x="X") +theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.x=element_text(face="italic"),plot.title = element_text(hjust = 0.5),axis.text = element_text(family="sans"),legend.text =element_text(size = 14)))

Eff_joint20=as.numeric(Samples_joint[,5])
eff_j=data.frame(Eff=Eff_joint20)
eff_j$Model<-'Joint'

eff_c=data.frame(Eff=Eff_cutposterior20)
eff_c$Model<-'Cut posterior'


count.df <- rbind(eff_j,eff_c)
plot <- ggplot(count.df, aes(Eff,fill=Model))
(pl2=plot + geom_histogram(aes(y=..density..),position = 'identity')+geom_vline(xintercept = -0.2/0.3,linetype='dotted',color='black') +xlim(min(append(Eff_joint20,Eff_cutposterior20)),max(append(Eff_joint20,Eff_cutposterior20)))+ 
  scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(title="Posterior probability of the treatment effect ",y="density",
                                                                              x="X") +theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))

ggarrange(pl1,pl2, common.legend = TRUE,legend = 'bottom')

## No treatment effect
# Sample the data

i=5000
set.seed(i)
X=sort(2*rbeta(500,2,2)-1)
Y=funC_sample(X)
T=sample_prob30(X)
{b_f1=bounds(X,25)
b_f1t=bounds(X,25)
ubr=b_f1$ubr
ubl=b_f1$ubl
ubrt=b_f1t$ubr
ublt=b_f1t$ubl
lb=b_f1$lb
b_s=bounds(X,50)
ubrs=b_s$ubr
ubls=b_s$ubl}
jlb=0.2
nc=1
# Set prior on c - uniform between -0.8 and 0.2
clb=-0.8
cub=0.2

dat1T=list(N=length(X),x=X,t=T,jlb=jlb,clb=-0.8,cub=0.2)
param_full=c('c','j','kl','kr','klt','krt','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t','tau1r','tau2r','tau1l','tau2l')
param_c=c('c')
initc1=list(c=-0.3,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
dat1_c<- run.jags('cutoff_initial_CONT.txt',inits = list(initc1) ,data=dat1T,monitor=param_c,burnin = 900,sample=2000,adapt = 100,n.chains = 1,method = 'simple')
C_start=as.numeric(combine.mcmc(dat1_c$mcmc))

#Fit treatment model
init1=Initial_treatment_CONT(X,T,C_start,lb,ubr,ubl,start,prob,100)
init2=Initial_treatment_CONT(X,T,C_start,lb,ubr,ubl,start,prob,200)
init3=Initial_treatment_CONT(X,T,C_start,lb,ubr,ubl,start,prob,300)
init4=Initial_treatment_CONT(X,T,C_start,lb,ubr,ubl,start,prob,400)

dat2=list(N=length(X),x=X,t=T,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,nc=nc,jlb=jlb,clb=clb,cub=cub)
dattreatment<-run.jags('treatment_CONT.txt', data=dat2,inits = list(init1,init2,init3,init4),monitor=param_full,burnin = 30000,sample=250,adapt = 500,n.chains = 4,method = 'parallel')

# Fit joint model
init1=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,100)  
init2=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,200)
init3=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,300)
init4=Initial_CONT_CONT(X,T,Y,C_start,lb,ubrs,ubls,jlb,400)
dat1=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubrt,ublt=ublt,lb=lb,nc=nc,jlb=jlb,clb=clb,cub=cub)
param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
system.time(datjoint<- run.jags('LoTTA_CONT_CONT', data=dat1,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 30000,sample=250,adapt = 500, method='parallel',n.chains = 4))

# Sample from cut posterior
Samples=combine.mcmc(dattreatment)
C=as.numeric(Samples[,1])
J=as.numeric(Samples[,2])
Eff_joint=as.numeric(Samples[,7])
hist(Eff_joint)
Eff_cutposterior=c()

for(i in 1:1000){
  print(i)
  c=C[i]
  init1=Initial_SHARP_CONT(X,Y,c,lb,ubrs,ubls,1)
  init2=Initial_SHARP_CONT(X,Y,c,lb,ubrs,ubls,2)
  init3=Initial_SHARP_CONT(X,Y,c,lb,ubrs,ubls,3)
  init4=Initial_SHARP_CONT(X,Y,c,lb,ubrs,ubls,4)
  dat2=list(N=length(X),x=X,y=Y,c=c,ubr=ubr,ubl=ubl,lb=lb)
  param=c('eff')
  datout<-run.jags('LoTTA_SHARP_CONT.txt', data=dat2,inits = list(init1,init2,init3,init4),monitor=param,burnin = 10000,sample=1,adapt = 500,n.chains = 1,method = 'parallel')
  Sampleseff=combine.mcmc(datout)
  eff=as.numeric(Sampleseff[1])/J[i]
  Eff_cutposterior=append(Eff_cutposterior,eff)
}

# Plot the results
Samples_joint=combine.mcmc(datjoint)
C_joint=as.numeric(Samples_joint[,1])
c_j=data.frame(C=C_joint)
c_j$Model<-'Joint'

c_c=data.frame(C=C)
c_c$Model<-'Cut posterior'

count.df <- rbind(c_j,c_c)
plot <- ggplot(count.df, aes(C,fill=Model))
(pl1=plot + geom_histogram(aes(y=..density..),position = 'identity') +xlim(-0.08,max(append(C_joint,C)))+ 
    scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5)))+geom_vline(xintercept = 0,linetype='dotted',color='black') +labs(title="Posterior probability of the cutoff location ",y="density",
                                                                                                                                           x="X") +theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.x=element_text(face="italic"),plot.title = element_text(hjust = 0.5),axis.text = element_text(family="sans"),legend.text =element_text(size = 14)))

Eff_joint=as.numeric(Samples_joint[,5])
eff_j=data.frame(Eff=Eff_joint)
eff_j$Model<-'Joint'

eff_c=data.frame(Eff=Eff_cutposterior)
eff_c$Model<-'Cut posterior'


count.df <- rbind(eff_j,eff_c)
plot <- ggplot(count.df, aes(Eff,fill=Model))
(pl2=plot + geom_histogram(aes(y=..density..),position = 'identity')+geom_vline(xintercept = 0,linetype='dotted',color='black') +xlim(min(append(Eff_joint,Eff_cutposterior)),max(append(Eff_joint,Eff_cutposterior)))+ 
    scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(title="Posterior probability of the treatment effect ",y="density",
                                                                                x="X") +theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))

ggarrange(pl1,pl2, common.legend = TRUE,legend = 'bottom')



