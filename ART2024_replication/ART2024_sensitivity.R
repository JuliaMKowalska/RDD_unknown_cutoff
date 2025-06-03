## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## ART Application - sensitivity analysis, Appendix B.4 ##
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
data <- read_dta("CKT_2023_SIM--ART.dta")
dd = data[complete.cases(data[,c("visit_test_6_18", "art_6m", "cd4")]),]
X=dd$cd4
Y=dd$visit_test_6_18
T=dd$art_6m
T=ifelse(T==1,0,1)

## Eta parameter

  # Prior on c - uniform between 300 and 399 (no change)
  start=299
  prob=rep(1,100) 
  # Lower bound on the jump size in the treatment probability function
  jlb=0.05 # eta 
  
  # Fitting two constant functions to initialize values of c (approx 50 seconds)
  datHIV_T=list(N=length(X),x=X,t=T,jlb=0.05,prob=prob,start=start,nc=nc,lb=lb)
  initcART=list(ct=50,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
  param_c=c('ct')
  datHIV_c<- run.jags('cutoff_initial_DIS.txt',inits = list(initcART) ,data=datHIV_T,monitor=param_c,burnin = 1000,sample=2000,adapt = 100,n.chains = 1,method = 'simple')
  Ct_start=as.numeric(combine.mcmc(datHIV_c))
  
  # Fitting treatment model (approx 35 minutes)
  init1=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,100,nc,jlb)
  init2=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,200,nc,jlb)
  init3=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,300,nc,jlb)
  init4=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,400,nc,jlb)
  param_cj=c('c','j')
  datHIV_T=list(N=length(X),x=X,t=T,jlb=jlb,prob=prob,start=start,nc=nc,lb=lb,ubrt=ubrt,ublt=ublt)
  system.time(HIV_FULL_treatment_eta<- run.jags('treatment_DIS.txt', data=datHIV_T,monitor=param_cj,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))

  
  # Fitting full LoTTA model (approx 13.5 hours)
  init1=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,100,nc,jlb)
  init2=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,200,nc,jlb)
  init3=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,300,nc,jlb)
  init4=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,400,nc,jlb)
  datHivfull=list(N=length(X),x=X,t=T,y=Y,ubr=b_hiv$ubr,ubl=b_hiv$ubl,ubrt=ubrt,ublt=ublt,lb=b_hiv$lb,prob=prob,start=start,nc=nc,jlb=jlb)
  param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
  system.time(HIV_FULL_eta<- run.jags('LoTTA_DIS_BIN', data=datHivfull,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))

## Informative prior
  
  # Prior on c - informative, p(c=350)=0.9
  start=299
  prob=rep(0.1,100) 
  prob[51]=0.1+0.9*100
  # Lower bound on the jump size in the treatment probability function (no change)
  jlb=0.2 # eta 
  
  # Fitting two constant functions to initialize values of c (approx 50 seconds)
  datHIV_T=list(N=length(X),x=X,t=T,jlb=0.05,prob=prob,start=start,nc=nc,lb=lb)
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
  system.time(HIV_FULL_treatment_inf<- run.jags('treatment_DIS.txt', data=datHIV_T,monitor=param_cj,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))
  
  # Fitting full LoTTA model (approx 13.5 hours)
  init1=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,100,nc)
  init2=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,200,nc)
  init3=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,300,nc)
  init4=Initial_DIS_BIN(X,T,Y,Ct_start,lb,ubr,ubl,start,prob,400,nc)
  datHivfull=list(N=length(X),x=X,t=T,y=Y,ubr=b_hiv$ubr,ubl=b_hiv$ubl,ubrt=ubrt,ublt=ublt,lb=b_hiv$lb,prob=prob,start=start,nc=nc,jlb=jlb)
  param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
  system.time(HIV_FULL_inf<- run.jags('LoTTA_DIS_BIN', data=datHivfull,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))


## Eta 0.05 vs Eta 0.2

# Load results from ART2024_main
load('HIV_FULL')
load('HIV_FULL_treatment')

HIV_FULL_treatment_05=HIV_FULL_treatment
#-------
Samples_FULL=combine.mcmc(HIV_FULL)
Samples_FULL_treatment=combine.mcmc(HIV_FULL_treatment)

Samples_FULL05=combine.mcmc(HIV_FULL_eta)
Samples_FULL_treatment05=combine.mcmc(HIV_FULL_treatment_eta)

# Figures

# -- Cutoff location --

  C_Full=as.numeric(Samples_FULL[,1])*nc
  C_Full05=as.numeric(Samples_FULL05[,1])*nc
  
  c_Full=data.frame(C=C_Full)  
  c_Full$Dataset<-'eta=0.2'
  
  c_Full05=data.frame(C=C_Full05)  
  c_Full05$Dataset<-'eta=0.05'
  
  count.df <- rbind(c_Full05,c_Full)
  plot <- ggplot(count.df, aes(C,fill=Dataset))
  (pl1=plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'dodge',binwidth = 0.5) +xlim(350,361)+ ylim(0,0.8)+
      scale_fill_manual(
        values = c(alpha("#E69F00", 0.6), alpha("black", 0.5)),
        labels = c(
          bquote(eta == 0.05),
          bquote(eta == 0.2)
        )) + labs(y="Density", x=expression(paste(italic("c "), "- cutoff" )),title = 'Posterior of the cutoff location, joint model',fill = "Lower bound") +
      theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))
  
  C_Treatment=as.numeric(Samples_FULL_treatment[,1])*nc
  C_Treatment05=as.numeric(Samples_FULL_treatment05[,1])*nc
  
  c_Treatment=data.frame(C=C_Treatment)  
  c_Treatment$Dataset<-'eta=0.2'
  
  c_Treatment05=data.frame(C=C_Treatment05)  
  c_Treatment05$Dataset<-'eta=0.05'
  
  count.df <- rbind(c_Treatment05,c_Treatment)
  plot <- ggplot(count.df, aes(C,fill=Dataset))
  (pl2=plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'dodge',binwidth = 0.5) +xlim(350,361)+ ylim(0,0.8)+
      scale_fill_manual(
        values = c(alpha("#E69F00", 0.6), alpha("black", 0.5)),
        labels = c(
          bquote(eta == 0.05),
          bquote(eta == 0.2)
        )) + labs(y="Density", x=expression(paste(italic("c "), "- cutoff" )),title = 'Posterior of the cutoff location, treatment model',fill = "Lower bound") +
      theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))
  plot=ggarrange(pl1,pl2,legend = 'bottom',common.legend = TRUE)
  ggsave("Eta_sensitivity_cutoff.pdf", plot = plot, width = 10, height = 5)

# -- Compliance rate --
  
  J_Full=as.numeric(Samples_FULL[,2])
  J_Full05=as.numeric(Samples_FULL05[,2])
  
  j_Full=data.frame(J=J_Full)  
  j_Full$Dataset<-'eta=0.2'
  
  j_Full05=data.frame(J=J_Full05)  
  j_Full05$Dataset<-'eta=0.05'
  ci(J_Treatment05,method='HDI',ci=0.95)
  
  count.df <- rbind(j_Full05,j_Full)
  plot <- ggplot(count.df, aes(J,fill=Dataset))
  (pl1=plot + geom_histogram(aes(y=..density..),position = 'identity') +xlim(min(append(J_Full,J_Full05)),max(append(J_Full,J_Full05)))+ 
      scale_fill_manual(
        values = c(alpha("#E69F00", 0.6), alpha("black", 0.5)),
        labels = c(
          bquote(eta == 0.05),
          bquote(eta == 0.2)
        )) +labs(
          y="Density", x=expression(paste(italic("j "), "- compliance rate" )),fill = "Prior",title = 'Posterior of the compliance rate, joint model') +
      theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))
  
  
  J_Treatment=as.numeric(Samples_FULL_treatment[,2])
  J_Treatment05=as.numeric(Samples_FULL_treatment05[,2])
  
  j_Treatment=data.frame(J=J_Treatment)  
  j_Treatment$Dataset<-'eta=0.2'
  
  j_Treatment05=data.frame(J=J_Treatment05)  
  j_Treatment05$Dataset<-'eta=0.05'
  
  
  count.df <- rbind(j_Treatment05,j_Treatment)
  plot <- ggplot(count.df, aes(J,fill=Dataset))
  (pl2=plot + geom_histogram(aes(y=..density..),position = 'identity') +xlim(min(append(J_Treatment,J_Treatment05)),max(append(J_Treatment,J_Treatment05)))+ 
      scale_fill_manual(
        values = c(alpha("#E69F00", 0.6), alpha("black", 0.5)),
        labels = c(
          bquote(eta == 0.05),
          bquote(eta == 0.2)
        )) +labs(
          y="Density", x=expression(paste(italic("j "), "- compliance rate" )),fill = "Prior",title = 'Posterior of the compliance rate, treatment model') +
      theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))
  
  plot=ggarrange(pl1,pl2,legend = 'bottom',common.legend = TRUE)
  ggsave("Eta_sensitivity_compliance.pdf", plot = plot, width = 10, height = 5)
  
# -- Treatment effect --
  
  Eff_Full=as.numeric(Samples_FULL[,5])
  Eff_Full05=as.numeric(Samples_FULL05[,5])
  
  # eta=0.2
  (eff_map_full=map_estimate(Eff_Full)) # MAP estimate
  (eff_hdi_full=ci(Eff_Full,method='HDI')) # HDI
  
  # eta=0.05
  (eff_map_full05=map_estimate(Eff_Full05)) # MAP estimate
  (eff_hdi_full05=ci(Eff_Full05,method='HDI')) # HDI
  
## Informative vs Uniform prior
  
Samples_FULL90=combine.mcmc(HIV_FULL_inf)
Samples_FULL_treatment90=combine.mcmc(HIV_FULL_treatment_inf) 

# Figures

# -- Cutoff location --

  C_Full=as.numeric(Samples_FULL[,1])*nc
  C_Full90=as.numeric(Samples_FULL90[,1])*nc
  
  c_Full=data.frame(C=C_Full)  
  c_Full$Dataset<-'Uniform'
  
  c_Full90=data.frame(C=C_Full90)  
  c_Full90$Dataset<-'Informative'
  
  count.df <- rbind(c_Full90,c_Full)
  plot <- ggplot(count.df, aes(C,fill=Dataset))
  (pl1=plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'dodge',binwidth = 0.5) +xlim(348,361)+ ylim(0,0.8)+
      scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) + labs(y="Density", x=expression(paste(italic("c "), "- cutoff" )),title = 'Posterior of the cutoff location, joint model',fill = "Prior") +
      theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))
  
  C_Treatment=as.numeric(Samples_FULL_treatment[,1])*nc
  C_Treatment90=as.numeric(Samples_FULL_treatment90[,1])*nc
  
  c_Treatment=data.frame(C=C_Treatment)  
  c_Treatment$Dataset<-'Uniform'
  
  c_Treatment90=data.frame(C=C_Treatment90)  
  c_Treatment90$Dataset<-'Informative'
  
  count.df <- rbind(c_Treatment90,c_Treatment)
  plot <- ggplot(count.df, aes(C,fill=Dataset))
  (pl2=plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'dodge',binwidth = 0.5) +xlim(348,361)+ ylim(0,0.8)+
      scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) + labs(y="Density", x=expression(paste(italic("c "), "- cutoff" )),title = 'Posterior of the cutoff location, treatment model',fill = "Lower bound") +
      theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))
  plot=ggarrange(pl1,pl2,legend = 'bottom',common.legend = TRUE)
  ggsave("C_sensitivity_cutoff.pdf", plot = plot, width = 10, height = 5)

# -- Compliance rate --
  
  J_Full=as.numeric(Samples_FULL[,2])
  J_Full90=as.numeric(Samples_FULL90[,2])
  
  j_Full=data.frame(J=J_Full)  
  j_Full$Dataset<-'Unifrom'
  
  j_Full90=data.frame(J=J_Full90)  
  j_Full90$Dataset<-'Informative'
  
  
  count.df <- rbind(j_Full90,j_Full)
  plot <- ggplot(count.df, aes(J,fill=Dataset))
  (pl1=plot + geom_histogram(aes(y=..density..),position = 'identity') +xlim(min(append(J_Full,J_Full90)),max(append(J_Full,J_Full90)))+ 
      scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
        y="Density", x=expression(paste(italic("j "), "- compliance rate" )),fill = "Prior",title = 'Posterior of the compliance rate, joint model') +
      theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))
  
  
  J_Treatment=as.numeric(Samples_FULL_treatment[,2])
  J_Treatment90=as.numeric(Samples_FULL_treatment90[,2])
  
  j_Treatment=data.frame(J=J_Treatment)  
  j_Treatment$Dataset<-'Uniform'
  
  j_Treatment90=data.frame(J=J_Treatment90)  
  j_Treatment90$Dataset<-'Informative'
  
  
  count.df <- rbind(j_Treatment90,j_Treatment)
  plot <- ggplot(count.df, aes(J,fill=Dataset))
  (pl2=plot + geom_histogram(aes(y=..density..),position = 'identity') +xlim(min(append(J_Treatment,J_Treatment90)),max(append(J_Treatment,J_Treatment90)))+ 
      scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
        y="Density", x=expression(paste(italic("j "), "- compliance rate" )),fill = "Prior",title = 'Posterior of the compliance rate, treatment model') +
      theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),plot.title = element_text(hjust = 0.5),legend.position = 'bottom'))
  
  plot=ggarrange(pl1,pl2,legend = 'bottom',common.legend = TRUE)
  ggsave("C_sensitivity_compliance.pdf", plot = plot, width = 10, height = 5)

# -- Treatment effect --
  
  Eff_Full=as.numeric(Samples_FULL[,5])
  Eff_Full90=as.numeric(Samples_FULL90[,5])
  
  # Uniform
  (eff_map_full=map_estimate(Eff_Full)) # MAP
  (eff_hdi_full=ci(Eff_Full,method='HDI')) # HDI
  
  # Informative
  (eff_map_full90=map_estimate(Eff_Full90)) # MAP
  (eff_hdi_full90=ci(Eff_Full90,method='HDI')) # HDI
  

