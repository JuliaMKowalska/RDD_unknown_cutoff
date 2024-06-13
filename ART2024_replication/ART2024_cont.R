## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## ART Application - additional code ##
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

data <- read.dta("CKT_2023_SIM--ART.dta")
dd = data[complete.cases(data[,c("visit_test_6_18", "art_6m", "cd4")]),]
X=dd$cd4
Y=dd$visit_test_6_18
T=dd$art_6m
T=ifelse(T==1,0,1)

## Main code ##
## The following code yields similar results to the code in ART2024_main.R that was used
# to obtain results in the paper, but it requires less computation time. Instead of
# a discrete prior on c, a continous prior is used. ## 
# Scatter plot of the treatment data
binsize=5
c=NULL
X_plot=Bin_data(T,X,c,binsize)[[2]]
T_plot=Bin_data(T,X,c,binsize)[[1]]
plot.df=data.frame('X'=X_plot,'Y'=T_plot)

ggscatter(plot.df, x = "X", y = "Y",color = "gray", size = 3, alpha = 0.6)+labs( y="Proportion T=1", x="CD4 Count",title = 'Prob. of delayed ART initiation') +
  theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 13),legend.position = 'bottom')

# Scatter plot of the outcome data
binsize=5
c=NULL
X_plot=Bin_data(Y,X,c,binsize)[[2]]
Y_plot=Bin_data(Y,X,c,binsize)[[1]]
plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
ggscatter(plot.df, x = "X", y = "Y",color = "gray", size = 3, alpha = 0.6)+labs( y="Proportion Y=1", x="CD4 Count",title = 'Prob. of retention in care') +
  theme_classic(base_size = 14)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),legend.position = 'bottom',plot.title = element_text(hjust = 0.5))

# -- FULL DATASET --
# Initial trimming & normalizing the data
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

# Prior on c - continous uniform between 300 and 400 (on th normalized dataset between 0.3 and 0.4)
clb=300/nc
cub=400/nc
# Lower bound on the jump size in the treatment probability function
jlb=0.2

# Fitting two constant functions to initialize values of c (approx 7 seconds)
datHIV_T=list(N=length(X),x=X,t=T,jlb=0.2,clb=clb,cub=cub,lb=lb)
initcART=list(c=350/nc,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
param_c=c('c')
system.time(datHIV_c<- run.jags('cutoff_initial_CONT.txt',inits = list(initcART) ,data=datHIV_T,monitor=param_c,burnin = 1000,sample=2000,adapt = 100,n.chains = 1,method = 'simple'))
C_start=as.numeric(combine.mcmc(datHIV_c))

# Fitting treatment model (approx 9 minutes)
init1=Initial_CONT_treatment(X,T,C_start,lb,ubr,ubl,start,prob,100)
init2=Initial_CONT_treatment(X,T,C_start,lb,ubr,ubl,start,prob,200)
init3=Initial_CONT_treatment(X,T,C_start,lb,ubr,ubl,start,prob,300)
init4=Initial_CONT_treatment(X,T,C_start,lb,ubr,ubl,start,prob,400)
param_cj=c('c','j')
datHIV_T=list(N=length(X),x=X,t=T,jlb=jlb,clb=clb,cub=cub,lb=lb,ubrt=ubrt,ublt=ublt)
system.time(HIV_FULL_treatment<- run.jags('treatment_CONT.txt', data=datHIV_T,monitor=param_cj,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))

# Convergence check
plot(HIV_FULL_treatment,c('trace'))

# Fitting full LoTTA model 
init1=Initial_CONT_BIN(X,T,Y,C_start,lb,ubr,ubl,start,prob,100)
init2=Initial_CONT_BIN(X,T,Y,C_start,lb,ubr,ubl,start,prob,200)
init3=Initial_CONT_BIN(X,T,Y,C_start,lb,ubr,ubl,start,prob,300)
init4=Initial_CONT_BIN(X,T,Y,C_start,lb,ubr,ubl,start,prob,400)
datHivfull=list(N=length(X),x=X,t=T,y=Y,ubr=b_hiv$ubr,ubl=b_hiv$ubl,ubrt=ubrt,ublt=ublt,lb=b_hiv$lb,clb=clb,cub=cub,nc=nc,jlb=jlb)
param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
system.time(HIV_FULL<- run.jags('LoTTA_BIN', data=datHivfull,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))
# Convergence check
plot(HIV_FULL,c('trace'),c('eff','c','j'))
gelman.diag(HIV_FULL)
# Plots
Samples_FULL=combine.mcmc(HIV_FULL)
Samples_FULL_treatment=combine.mcmc(HIV_FULL_treatment)
# -- Cutoff location --
C_Full=as.numeric(Samples_FULL[,1])*nc
C_Full_DIS=ceiling(C_Full) # discritized cutoff posterior
C_Treatment=as.numeric(Samples_FULL_treatment[,1])*nc
C_Treatment_DIS=ceiling(C_Treatment) # discritized cutoff posterior

c_Full=data.frame(C=C_Full_DIS)  
c_Full$Model<-'Full'

c_Treatment=data.frame(C=C_Treatment_DIS)
c_Treatment$Model<-'Treatment'

count.df <- rbind(c_Treatment,c_Full)
plot <- ggplot(count.df, aes(C,fill=Model))
plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'dodge',binwidth = 0.5) +xlim(350,361)+ 
  scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
    y="Density", x=expression(paste(italic("c "), "- cutoff" )),title = 'Posterior of the cutoff location') +
  theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),legend.position = 'bottom')

# -- Compliance rate --
J_Full=as.numeric(Samples_FULL[,2])
J_Treatment=as.numeric(Samples_FULL_treatment[,2])

j_Full=data.frame(J=J_Full)  
j_Full$Model<-'Full'

j_Treatment=data.frame(J=J_Treatment)
j_Treatment$Model<-'Treatment'

count.df <- rbind(j_Treatment,j_Full)
plot <- ggplot(count.df, aes(J,fill=Model))
plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'identity') +xlim(min(append(J_Full,J_Treatment)),max(append(J_Full,J_Treatment)))+ 
  scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
    y="Density", x=expression(paste(italic("j "), "- compliance rate" )),title = 'Posterior of the compliance rate') +
  theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),legend.position = 'bottom')


# -- TRIMMED DATASET --
ddtr=subset(dd,dd$cd4>=200 & dd$cd4<=500)
X=ddtr$cd4/nc
Y=ddtr$visit_test_6_18
T=ddtr$art_6m
T=ifelse(T==1,0,1)

{b_hiv=bounds(X,25)
  clb=300
  cub=400
  ubr=b_hiv$ubr
  ubl=b_hiv$ubl
  lb=b_hiv$lb
  b_hivt=bounds(X,25)
  ubrt=b_hivt$ubr
  ublt=b_hivt$ubl}


datHIV_T=list(N=length(X),x=X,t=T,jlb=jlb,clb=clb,cub=cub,lb=lb)
initcART=list(ct=0.35,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
param_c=c('c')
datHIV_cTr<- run.jags('cutoff_initial_CONT.txt',inits = list(initcART) ,data=datHIV_T,monitor=param_c,burnin = 1000,sample=2000,adapt = 100,n.chains = 1,method = 'simple')
C_startTr=as.numeric(combine.mcmc(datHIV_cTr))

# Fitting full LoTTA model (approx 1.5 hours)
init1=Initial_CONT_BIN(X,T,Y,C_startTr,lb,ubr,ubl,start,prob,100)
init2=Initial_CONT_BIN(X,T,Y,C_startTr,lb,ubr,ubl,start,prob,200)
init3=Initial_CONT_BIN(X,T,Y,C_startTr,lb,ubr,ubl,start,prob,300)
init4=Initial_CONT_BIN(X,T,Y,C_startTr,lb,ubr,ubl,start,prob,400)

datHivTr=list(N=length(X),x=X,t=T,y=Y,ubr=b_hiv$ubr,ubl=b_hiv$ubl,ubrt=ubrt,ublt=ublt,lb=b_hiv$lb,clb=clb,cub=cub,nc=nc,jlb=jlb)
param_full=c('c','j','kl','kr','eff','a0l','a1l','a2l','a3l','a0r','a1r','a2r','a3r','b1lt','a1lt','a2lt','b2lt','b1rt','a1rt','a2rt','b2rt','k1t','k2t')
system.time(HIV_TRIM<- run.jags('LoTTA_CONT_BIN', data=datHivTr,monitor=param_full,inits =list(init1,init2,init3,init4),burnin = 40000,sample=25000,adapt = 1000, method='parallel',n.chains = 4))

# Convergence check
plot(HIV_TRIM,c('trace'),c('eff','c','j'))
gelman.diag(HIV_TRIM)
# Plots
Samples_TRIM=combine.mcmc(HIV_TRIM)
# -- Cutoff location --
C_Full=as.numeric(Samples_FULL[,1])*nc
C_Full_DIS=ceiling(C_Full)
C_Trimmed=as.numeric(Samples_TRIM[,1])*nc
C_Trimmed_DIS=ceiling(C_Trimmed)

c_Full=data.frame(C=C_Full_DIS)  
c_Full$Dataset<-'Full'

c_Trimmed=data.frame(C=C_Trimmed_DIS)
c_Trimmed$Dataset<-'Trimmed'

count.df <- rbind(c_Trimmed,c_Full)
plot <- ggplot(count.df, aes(C,fill=Dataset))
plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'dodge',binwidth = 0.5) +xlim(350,361)+ 
  scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
    y="Density", x=expression(paste(italic("c "), "- cutoff" )),title = 'Posterior of the cutoff location') +
  theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),legend.position = 'bottom')

# -- Compliance rate --
J_Full=as.numeric(Samples_FULL[,2])
J_Trimmed=as.numeric(Samples_TRIM[,2])

j_Full=data.frame(J=J_Full)  
j_Full$Dataset<-'Full'

j_Trimmed=data.frame(J=J_Trimmed)
j_Trimmed$Dataset<-'Trimmed'

count.df <- rbind(j_Trimmed,j_Full)
plot <- ggplot(count.df, aes(J,fill=Dataset))
plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'identity') +xlim(min(append(J_Full,J_Trimmed)),max(append(J_Full,J_Trimmed)))+ 
  scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
    y="Density", x=expression(paste(italic("j "), "- compliance rate" )),title = 'Posterior of the compliance rate') +
  theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),legend.position = 'bottom')

# -- Treatment effect --
Eff_Full=as.numeric(Samples_FULL[,5])
Eff_Trimmed=as.numeric(Samples_TRIM[,5])

eff_Full=data.frame(Eff=Eff_Full)  
eff_Full$Dataset<-'Full'

eff_Trimmed=data.frame(Eff=Eff_Trimmed)
eff_Trimmed$Dataset<-'Trimmed'

count.df <- rbind(eff_Trimmed,eff_Full)
plot <- ggplot(count.df, aes(Eff,fill=Dataset))
plot + geom_histogram(aes(y=2*after_stat(count)/sum(after_stat(count))),position = 'identity') +xlim(min(append(Eff_Full,Eff_Trimmed)),max(append(Eff_Full,Eff_Trimmed)))+ 
  scale_fill_manual(values=c(alpha("#E69F00",0.6), alpha("black",0.5))) +labs(
    y="Density", x=expression(paste(italic("eff "), "- treatment effect" )),title = 'Posterior of the treatment effect') +
  theme_classic(base_size = 12)+theme(text = element_text(family='serif'),axis.text=element_text(family = "sans"),legend.text =element_text(size = 14),legend.position = 'bottom')
# -- Posterior functions --
# Focus region - treatment probablity
lb_p=200
ub_p=500
ddtr=subset(dd,dd$cd4>=lb_p & dd$cd4<=ub_p)
X=ddtr$cd4
Y=ddtr$visit_test_6_18
T=ddtr$art_6m
T=ifelse(T==1,0,1)
binsize=10
c=355
X_plot=Bin_data(T,X,c,binsize)[[2]]
T_plot=Bin_data(T,X,c,binsize)[[1]]
x=seq(lb_p,ub_p,1)

fun_sampleF=apply(Samples_FULL,1,treatment_function_sample,x=x,nc=nc)
q_funF=apply(fun_sampleF,1,quantile,probs=c(0.025,0.5,0.975))

fun_sampleT=apply(Samples_TRIM,1,treatment_function_sample,x=x,nc=nc)
q_funT=apply(fun_sampleT,1,quantile,probs=c(0.025,0.5,0.975))

plot.df=data.frame('X'=X_plot,'Y'=T_plot)
lineF.df=data.frame(X=x,Y=q_funF[2,],lower95=q_funF[1,],upper95=q_funF[3,])
lineT.df=data.frame(X=x,Y=q_funT[2,],lower95=q_funT[1,],upper95=q_funT[3,])


ggscatter(plot.df, x = "X", y = "Y",
          color = "gray",
          size = 3, alpha = 0.6)+geom_vline(xintercept = 355,linetype='dotted',color="#009E73")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Full'),data=lineF.df,alpha=0.35)+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1,col='#E69F00')+geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Trimmed'),data=lineT.df,alpha=0.2)+geom_point(aes(x=X,y=Y),data=lineT.df,size=0.1)+labs( y="", x=expression(paste(italic('X')," - score")),title="Treatment probability function",subtitle = 'Prob. of delayed ART initiation')+theme_classic(base_size = 13)+theme(text = element_text(family='serif'),legend.text =element_text(size = 14),legend.position = 'bottom')+ scale_fill_manual(name="Model fit",values=c('Full'='#E69F00','Trimmed'='black'),labels = c(expression(paste("Full dataset ", italic(X) %in%"[50,950]"  )), expression(paste("Trimmed dataset ", italic(X) %in%"[200,500]"  )))) 

# Focus region - outcome function
X_plot=Bin_data(Y,X,c,binsize)[[2]]
Y_plot=Bin_data(Y,X,c,binsize)[[1]]

fun_sampleF=apply(Samples_FULL,1,binary_outcome_function_sample,x=x,nc=nc)
q_funF=apply(fun_sampleF,1,quantile,probs=c(0.025,0.5,0.975))

fun_sampleT=apply(Samples_TRIM,1,binary_outcome_function_sample,x=x,nc=nc)
q_funT=apply(fun_sampleT,1,quantile,probs=c(0.025,0.5,0.975))

plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=q_funF[2,],lower95=q_funF[1,],upper95=q_funF[3,])
lineT.df=data.frame(X=x,Y=q_funT[2,],lower95=q_funT[1,],upper95=q_funT[3,])

ggscatter(plot.df, x = "X", y = "Y",
          color = "gray",
          size = 3, alpha = 0.6)+geom_vline(xintercept = 355,linetype='dotted',color="#009E73")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Full'),data=lineF.df,alpha=0.35)+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1,col='#E69F00')+geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Trimmed'),data=lineT.df,alpha=0.2)+geom_point(aes(x=X,y=Y),data=lineT.df,size=0.1)+labs( y="", x=expression(paste(italic('X')," - score")),title = 'Outcome function',subtitle = 'Prob. of retention in care')+theme_classic(base_size = 13)+theme(text = element_text(family='serif'),legend.text =element_text(size = 14),legend.position = 'bottom')+ scale_fill_manual(name="Model fit",values=c('Full'='#E69F00','Trimmed'='black'),labels = c(expression(paste("Full dataset ", italic(X) %in%"[50,950]"  )), expression(paste("Trimmed dataset ", italic(X) %in%"[200,500]"  )))) 


# Full dataset - treatment probablity
lb_p=50
ub_p=950
ddtr=subset(dd,dd$cd4>=lb_p & dd$cd4<=ub_p)
X=ddtr$cd4
Y=ddtr$visit_test_6_18
T=ddtr$art_6m
T=ifelse(T==1,0,1)
binsize=10
c=355 
X_plot=Bin_data(T,X,c,binsize)[[2]]
T_plot=Bin_data(T,X,c,binsize)[[1]]
x=seq(lb_p,ub_p,1)

fun_sampleF=apply(Samples_FULL,1,treatment_function_sample,x=x,nc=nc)
q_funF=apply(fun_sampleF,1,quantile,probs=c(0.025,0.5,0.975))

fun_sampleT=apply(Samples_TRIM,1,treatment_function_sample,x=x,nc=nc)
q_funT=apply(fun_sampleT,1,quantile,probs=c(0.025,0.5,0.975))

plot.df=data.frame('X'=X_plot,'Y'=T_plot)
lineF.df=data.frame(X=x,Y=q_funF[2,],lower95=q_funF[1,],upper95=q_funF[3,])
lineT.df=data.frame(X=x,Y=q_funT[2,],lower95=q_funT[1,],upper95=q_funT[3,])


ggscatter(plot.df, x = "X", y = "Y",
          color = "gray",
          size = 3, alpha = 0.6)+geom_vline(xintercept = 355,linetype='dotted',color="#009E73")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Full'),data=lineF.df,alpha=0.35)+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1,col='#E69F00')+labs( y="", x=expression(paste(italic('X')," - score")),title="Treatment probability function",subtitle = 'Prob. of delayed ART initiation')+theme_classic(base_size = 13)+theme(text = element_text(family='serif'),legend.text =element_text(size = 14),legend.position = 'bottom')+ scale_fill_manual(name="Model fit",values=c('Full'='#E69F00','Trimmed'='black'),labels = c(expression(paste("Full dataset ", italic(X) %in%"[50,950]"  )), expression(paste("Trimmed dataset ", italic(X) %in%"[200,500]"  )))) 

# Focus region - outcome function
X_plot=Bin_data(Y,X,c,binsize)[[2]]
Y_plot=Bin_data(Y,X,c,binsize)[[1]]

fun_sampleF=apply(Samples_FULL,1,binary_outcome_function_sample,x=x,nc=nc)
q_funF=apply(fun_sampleF,1,quantile,probs=c(0.025,0.5,0.975))

fun_sampleT=apply(Samples_TRIM,1,binary_outcome_function_sample,x=x,nc=nc)
q_funT=apply(fun_sampleT,1,quantile,probs=c(0.025,0.5,0.975))

plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=q_funF[2,],lower95=q_funF[1,],upper95=q_funF[3,])
lineT.df=data.frame(X=x,Y=q_funT[2,],lower95=q_funT[1,],upper95=q_funT[3,])

ggscatter(plot.df, x = "X", y = "Y",
          color = "gray",
          size = 3, alpha = 0.6)+geom_vline(xintercept = 355,linetype='dotted',color="#009E73")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_ribbon(aes(x=X,y=Y,ymin = lower95, ymax = upper95,fill='Full'),data=lineF.df,alpha=0.35)+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1,col='#E69F00')+labs( y="", x=expression(paste(italic('X')," - score")),title = 'Outcome function',subtitle = 'Prob. of retention in care')+theme_classic(base_size = 13)+theme(text = element_text(family='serif'),legend.text =element_text(size = 14),legend.position = 'bottom')+ scale_fill_manual(name="Model fit",values=c('Full'='#E69F00','Trimmed'='black'),labels = c(expression(paste("Full dataset ", italic(X) %in%"[50,950]"  )), expression(paste("Trimmed dataset ", italic(X) %in%"[200,500]"  )))) 

# LLR approach
# Full dataset c=355
ddtr=subset(dd,dd$cd4>=50 & dd$cd4<=950)
X=ddtr$cd4
Y=ddtr$visit_test_6_18
T=ddtr$art_6m
T=ifelse(T==1,0,1)
r_FULL=rdrobust(Y,X,355,fuzzy=T)
eff_LLR_full=r_FULL$Estimate[2]
eff_ci_LLR_full=as.numeric(r_FULL$ci[3,])
j_LLR_full=r_FULL$tau_T[2]
j_ci_LLR_full=as.numeric(r_FULL$ci_T[3,])

# Trimmed dataset c=355
ddtr=subset(dd,dd$cd4>=200 & dd$cd4<=500)
X=ddtr$cd4
Y=ddtr$visit_test_6_18
T=ddtr$art_6m
T=ifelse(T==1,0,1)
r_TRIM=rdrobust(Y,X,355,fuzzy=T)
eff_LLR_trim=r_TRIM$Estimate[2]
eff_ci_LLR_trim=as.numeric(r_TRIM$ci[3,])
j_LLR_trim=r_TRIM$tau_T[2]
j_ci_LLR_trim=as.numeric(r_TRIM$ci_T[3,])

# Full dataset c=350
ddtr=subset(dd,dd$cd4>=50 & dd$cd4<=950)
X=ddtr$cd4
Y=ddtr$visit_test_6_18
T=ddtr$art_6m
T=ifelse(T==1,0,1)
r_FULL=rdrobust(Y,X,350,fuzzy=T)
eff_LLR_full2=r_FULL$Estimate[2]
eff_ci_LLR_full2=as.numeric(r_FULL$ci[3,])
j_LLR_full2=r_FULL$tau_T[2]
j_ci_LLR_full2=as.numeric(r_FULL$ci_T[3,])

# Trimmed dataset c=350
ddtr=subset(dd,dd$cd4>=200 & dd$cd4<=500)
X=ddtr$cd4
Y=ddtr$visit_test_6_18
T=ddtr$art_6m
T=ifelse(T==1,0,1)
r_TRIM=rdrobust(Y,X,350,fuzzy=T)
eff_LLR_trim2=r_TRIM$Estimate[2]
eff_ci_LLR_trim2=as.numeric(r_TRIM$ci[3,])
j_LLR_trim2=r_TRIM$tau_T[2]
j_ci_LLR_trim2=as.numeric(r_TRIM$ci_T[3,])

# LoTTA Estimates

# Treatment effect
eff_map_full=map_estimate(Eff_Full)
eff_hdi_full=ci(Eff_Full,method='HDI')
eff_map_trimmed=map_estimate(Eff_Trimmed)
eff_hdi_trimmed=ci(Eff_Trimmed,method='HDI')
paste("FULL DATASET: MAP estimate of the treatment effect with 95% hdi interval:", round(eff_map_full,2),"(",round(as.numeric(eff_hdi_full[2]),2),",",round(as.numeric(eff_hdi_full[3]),2),")")
paste("FULL DATASET: LLR estimate of the treatment effect at c=355 with 95% conf. interval:", round(eff_LLR_full,2),"(",round(as.numeric(eff_ci_LLR_full[1]),2),",",round(as.numeric(eff_ci_LLR_full[2]),2),")")
paste("FULL DATASET: LLR estimate of the treatment effect at c=350 with 95% conf. interval:", round(eff_LLR_full2,2),"(",round(as.numeric(eff_ci_LLR_full2[1]),2),",",round(as.numeric(eff_ci_LLR_full2[2]),2),")")

paste("TRIMMED DATASET: MAP estimate of the treatment effect with 95% hdi interval:", round(eff_map_trimmed,2),"(",round(as.numeric(eff_hdi_trimmed[2]),2),",",round(as.numeric(eff_hdi_trimmed[3]),2),")")
paste("TRIMMED DATASET: LLR estimate of the treatment effect at c=355 with 95% conf. interval:", round(eff_LLR_trim,2),"(",round(as.numeric(eff_ci_LLR_trim[1]),2),",",round(as.numeric(eff_ci_LLR_trim[2]),2),")")
paste("TRIMMED DATASET: LLR estimate of the treatment effect at c=350 with 95% conf. interval:", round(eff_LLR_trim2,2),"(",round(as.numeric(eff_ci_LLR_trim2[1]),2),",",round(as.numeric(eff_ci_LLR_trim2[2]),2),")")

# Complinace rate

j_map_full=map_estimate(J_Full)
j_hdi_full=ci(J_Full,method='HDI')
j_map_trimmed=map_estimate(J_Trimmed)
j_hdi_trimmed=ci(J_Trimmed,method='HDI')
paste("FULL DATASET: MAP estimate of the compliance rate with 95% hdi interval:", round(j_map_full,2),"(",round(as.numeric(j_hdi_full[2]),2),",",round(as.numeric(j_hdi_full[3]),2),")")
paste("FULL DATASET: LLR estimate of the compliance rate at c=355 with 95% conf. interval:", round(j_LLR_full,2),"(",round(as.numeric(j_ci_LLR_full[1]),2),",",round(as.numeric(j_ci_LLR_full[2]),2),")")
paste("FULL DATASET: LLR estimate of the compliance rate at c=350 with 95% conf. interval:", round(j_LLR_full2,2),"(",round(as.numeric(j_ci_LLR_full2[1]),2),",",round(as.numeric(j_ci_LLR_full2[2]),2),")")

paste("TRIMMED DATASET: MAP estimate of the compliance rate with 95% hdi interval:", round(j_map_trimmed,2),"(",round(as.numeric(j_hdi_trimmed[2]),2),",",round(as.numeric(j_hdi_trimmed[3]),2),")")
paste("TRIMMED DATASET: LLR estimate of the compliance rate at c=355 with 95% conf. interval:", round(j_LLR_trim,2),"(",round(as.numeric(j_ci_LLR_trim[1]),2),",",round(as.numeric(j_ci_LLR_trim[2]),2),")")
paste("TRIMMED DATASET: LLR estimate of the compliance rate at c=350 with 95% conf. interval:", round(j_LLR_trim2,2),"(",round(as.numeric(j_ci_LLR_trim2[1]),2),",",round(as.numeric(j_ci_LLR_trim2[2]),2),")")

# Cutoff

c_map_full=map_estimate(C_Full)
c_map_trimmed=map_estimate(C_Trimmed)
paste("FULL DATASET: MAP estimate of the cutoff:", round(c_map_full,2))
paste("TRIMMED DATASET: MAP estimate of the cutoff:", round(c_map_trimmed,2))
