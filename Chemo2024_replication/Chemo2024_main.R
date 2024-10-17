## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## Chemotherapy Application - main code ##
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

source('Chemo2024_functions.R')

## The data was made available by Matias D. Cattaneo and can be downloaded under the link:
# https://github.com/rdpackages-replication/CKT_2023_SIM/blob/master/CKT_2023_Chemo.dta ##

#download.file("https://github.com/rdpackages-replication/CKT_2023_SIM/blob/master/CKT_2023_Chemo.dta",destfile = "CKT_2023_Chemo.dta")
data <- read.dta("CKT_2023_Chemo.dta", warn.missing.labels = FALSE)
X=data$onc_score
T=data$chemo
Y=data$cancer2

# Scatter plot of the treatment data
binsize=2
c=NULL
X_plot=Bin_data(T,X,c,binsize)[[2]]
T_plot=Bin_data(T,X,c,binsize)[[1]]

plot.df=data.frame(X=X_plot,Y=T_plot)
x.title<-expression(paste(italic('X')," - score, ", italic("c")," = 0 ", "cutoff"))
(plT=ggscatter(plot.df, x = "X", y = "Y",
               color = "gray",
               size = 5, alpha = 0.6)+labs(title="Prob. of receiving chemotherapy",
                                           x ='Oncotype score', y = expression(paste("Proportion ","T=1")))+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),plot.title = element_text(hjust = 0.5),axis.text=element_text(family = "sans")))

# Scatter plot of the outcome data
Y_plot=Bin_data(Y,X,c,binsize)[[1]]
plot.df=data.frame(X=X_plot,Y=Y_plot)
(plY=ggscatter(plot.df, x = "X", y = "Y",
               color = "gray",
               size = 5, alpha = 0.6)+labs(title="Prob. of cancer recurrence",
                                           x ='Oncotype score', y = expression(paste("Proportion ","Y=1")))+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),plot.title = element_text(hjust = 0.5),axis.text=element_text(family = "sans")))
ggarrange(plT,plY)
# Normalizing the data
nc=100 # normalizing constant
X=data$onc_score/nc

# Window bounds
{b_f1=bounds(X,25)
  ubr=b_f1$ubr
  ubl=b_f1$ubl
  lb=b_f1$lb
  b_s=bounds(X,50)
  ubrs=b_s$ubr
  ubls=b_s$ubl}

# Prior on c - discrete uniform between 20 and 30
start=19
prob=rep(1,11)

# Fit treatment and full LoTTA model with jlb=0.2
jlb=0.2
# Fit constant model
datchemo_T=list(N=length(X),x=X,t=T,jlb=jlb,prob=prob,start=start,nc=nc)
initcCHEMO=list(ct=10,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
param_c=c('ct')
datchemo_c<- run.jags('cutoff_initial_DIS.txt' ,data=datchemo_T,monitor=param_c,inits=list(initcCHEMO),burnin = 1000,sample=1000,adapt = 1000,n.chains = 1,method = 'parallel')
Ct_start=as.numeric(combine.mcmc(datchemo_c))
# Fit treatment model (1.2 min)
init1=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,1,nc)
init2=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,2,nc)
init3=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,3,nc)
init4=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,4,nc)
param_cj=c('c','j')
dat1=list(N=length(X),x=X,t=T,jlb=jlb,prob=prob,start=start,nc=nc,lb=lb,ubrt=ubr,ublt=ubl)
system.time(Chemo_treatment20<- run.jags('treatment_DIS.txt', data=dat1,monitor=param_cj,inits =list(init1,init2,init3,init4),burnin = 30000,sample=25000,adapt = 500, method='parallel',n.chains = 4))
# Fit full model (40 min)
init1=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,1,nc)
init2=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,2,nc)
init3=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,3,nc)
init4=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,4,nc)

dat2=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubr,ublt=ubl,lb=lb,prob=prob,start=start,nc=nc,jlb=jlb)
param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')

system.time(Chemo_final20<- run.jags('LoTTA_DIS_BIN.txt', data=dat2,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 30000,sample=25000,adapt = 500, method='parallel',n.chains = 4))

# Fit treatment and full LoTTA model with jlb=0.1
jlb=0.1
# Fit constant model
datchemo_T=list(N=length(X),x=X,t=T,jlb=jlb,prob=prob,start=start,nc=nc)
initcCHEMO=list(ct=10,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
param_c=c('ct')
datchemo_c<- run.jags('cutoff_initial_DIS.txt' ,data=datchemo_T,monitor=param_c,inits=list(initcCHEMO),burnin = 1000,sample=1000,adapt = 1000,n.chains = 1,method = 'parallel')
Ct_start=as.numeric(combine.mcmc(datchemo_c))
# Fit treatment model (1.2 min)
init1=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,1,nc,jlb)
init2=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,2,nc,jlb)
init3=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,3,nc,jlb)
init4=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,4,nc,jlb)
param_cj=c('c','j')
dat1=list(N=length(X),x=X,t=T,jlb=jlb,prob=prob,start=start,nc=nc,lb=lb,ubrt=ubr,ublt=ubl)
system.time(Chemo_treatment10<- run.jags('treatment_DIS.txt', data=dat1,monitor=param_cj,inits =list(init1,init2,init3,init4),burnin = 30000,sample=25000,adapt = 500, method='parallel',n.chains = 4))
# Fit full model (40 minutes)
init1=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,1,nc,jlb)
init2=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,2,nc,jlb)
init3=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,3,nc,jlb)
init4=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,4,nc,jlb)

dat2=list(N=length(X),x=X,t=T,y=Y,ubr=ubr,ubl=ubl,ubrt=ubr,ublt=ubl,lb=lb,prob=prob,start=start,nc=nc,jlb=jlb)
param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')

system.time(Chemo_final10<- run.jags('LoTTA_DIS_BIN.txt', data=dat2,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 30000,sample=25000,adapt = 500, method='parallel',n.chains = 4))

# Plots
Samples_FULL20=combine.mcmc(Chemo_final20)
Samples_FULL10=combine.mcmc(Chemo_final10)
Samples_TREATMENT20=combine.mcmc(Chemo_treatment20)
Samples_TREATMENT10=combine.mcmc(Chemo_treatment10)


C_Full20=as.numeric(Samples_FULL20[,1])*nc
C_Treatment20=as.numeric(Samples_TREATMENT20[,1])*nc

c_Full20=data.frame(C=C_Full20)  
c_Full20$Model<-'Full'

c_Treatment20=data.frame(C=C_Treatment20)
c_Treatment20$Model<-'Treatment'

count.df <- rbind(c_Treatment20,c_Full20)
plot <- ggplot(count.df, aes(C,fill=Model))
(pl1=plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'dodge',binwidth = 0.5) +xlim(23,30)+ 
    scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
      y="Density", x=expression(paste(italic("c "), "- cutoff" )),title = 'Posterior of the cutoff location (j>0.2)') +
    theme_classic(base_size = 14)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))


J_Full20=as.numeric(Samples_FULL20[,2])
J_Treatment20=as.numeric(Samples_TREATMENT20[,2])

j_Full20=data.frame(J=J_Full20)  
j_Full20$Model<-'Full'

j_Treatment20=data.frame(J=J_Treatment20)
j_Treatment20$Model<-'Treatment'

count.df <- rbind(j_Treatment20,j_Full20)
plot <- ggplot(count.df, aes(J,fill=Model))
(pl2=plot + geom_histogram(aes(y=..density..),position = 'identity') +xlim(min(append(J_Full20,J_Treatment20)),max(append(J_Full20,J_Treatment20)))+ 
    scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
      y="Density", x=expression(paste(italic("j "), "- compliance rate" )),title = 'Posterior of the compliance rate (j>0.2)') +
    theme_classic(base_size = 14)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))

ggarrange(pl1,pl2,common.legend = TRUE,legend="bottom")

C_Full10=as.numeric(Samples_FULL10[,1])*nc
C_Treatment10=as.numeric(Samples_TREATMENT10[,1])*nc

c_Full10=data.frame(C=C_Full10)  
c_Full10$Model<-'Full'

c_Treatment10=data.frame(C=C_Treatment10)
c_Treatment10$Model<-'Treatment'

count.df <- rbind(c_Treatment10,c_Full10)
plot <- ggplot(count.df, aes(C,fill=Model))
(pl3=plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'dodge',binwidth = 0.5) +xlim(23,30)+ 
    scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
      y="Density", x=expression(paste(italic("c "), "- cutoff" )),title = 'Posterior of the cutoff location (j>0.1)') +
    theme_classic(base_size = 14)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))


J_Full10=as.numeric(Samples_FULL10[,2])
J_Treatment10=as.numeric(Samples_TREATMENT10[,2])

j_Full10=data.frame(J=J_Full10)  
j_Full10$Model<-'Full'

j_Treatment10=data.frame(J=J_Treatment10)
j_Treatment10$Model<-'Treatment'

count.df <- rbind(j_Treatment10,j_Full10)
plot <- ggplot(count.df, aes(J,fill=Model))
(pl4=plot + geom_histogram(aes(y=..density..),position = 'identity') +xlim(min(append(J_Full10,J_Treatment10)),max(append(J_Full20,J_Treatment20)))+ 
    scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
      y="Density", x=expression(paste(italic("j "), "- compliance rate" )),title = 'Posterior of the compliance rate (j>0.1)') +
    theme_classic(base_size = 14)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))

ggarrange(pl1,pl2,pl3,pl4,common.legend = TRUE,legend="bottom")
