## Replication code for "Bayesian Regression Discontinuity Design with Unknown Cutoff" ##
## Simulations - main script & supplementary materials ##
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

source('Simulations2024_functions.R')

# The following statistics are calculated:
# LoTTA:
# abs_err_med, abs_err_map - absolute error of median/MAP estimate of the treatment effect, 
# ci_length_sym, ci_length_hdi - length of symmetric/HDI 95% credible interval
# bias_med, bias_map - bias of median/MAP estimate of the treatment effect, 
# cov_med, cov_hdi - coverage of symmetric/HDI 95% credible interval,
# sign_med, sign_hdi - correct sign identification of symmetric/HDI 95% credible interval,
# c_abs_med, c_abs_map - absolute error of median/MAP estimate of the cutoff, 
# j_abs_med, j_abs_map - absolute error of median/MAP estimate of the compliance rate
# LLR:
# abs_err - absolute error of bias corrected estimate of the treatment effect, 
# ci_length - length of robust 95% confidence interval
# bias - bias of bias corrected estimate of the treatment effect, 
# cov - coverage of robust 95% confidence interval,
# sign - correct sign identification of robust 95% confidence interval,
# j_abs - absolute error of the estimate of the compliance rate

#----------------------------#
#  Section 5 & Appendix B.1  #
#----------------------------#

## Fuzzy design ##
# A single simulation takes around 2.2 minutes. In the following code we split simulations
# into two batches of 500.

# Function A #

# jump = 0.55

system.time(postA_1to500<-lapply(2001:2500, simulation_FUZZY,"A","0.55"))
system.time(postA_501to1000<-lapply(2501:3000, simulation_FUZZY,"A","0.55"))
# Results
  Results1<-lapply(postA_1to500, performance_sample_FUZZY,"A","0.55")
  Results2<-lapply(postA_501to1000, performance_sample_FUZZY,"A","0.55")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results[,1:6], 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)
  sqrt(mean(Results[['c_abs_map']]^2)) # RMSE cutoff location (MAP estimate)
  sqrt(mean(Results[['j_abs_map']]^2)) # RMSE compliance rate (MAP estimate)




# LLR
LLR_res=LLR_performance_FUZZY(1000,'A','0.55',2000)
apply(LLR_res, 2,mean)
apply(LLR_res[1:3], 2,median)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. 
sqrt(mean(LLR_res[['j_abs']]^2)) # RMSE compliance rate 

 
# jump = 0.3

system.time(postA_1to500<-lapply(2001:2500, simulation_FUZZY,"A","0.3"))
system.time(postA_501to1000<-lapply(2501:3000, simulation_FUZZY,"A","0.3"))

# Results
  Results1<-lapply(postA_1to500, performance_sample_FUZZY,"A","0.3")
  Results2<-lapply(postA_501to1000, performance_sample_FUZZY,"A","0.3")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results[,1:6], 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)
  sqrt(mean(Results[['c_abs_map']]^2)) # RMSE cutoff location (MAP estimate)
  sqrt(mean(Results[['j_abs_map']]^2)) # RMSE compliance rate (MAP estimate)


# LLR
LLR_res=LLR_performance_FUZZY(1000,'A','0.3',2000)
apply(LLR_res, 2,mean)
apply(LLR_res[1:3], 2,median)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. 
sqrt(mean(LLR_res[['j_abs']]^2)) # RMSE compliance rate 


# Function B #

# jump = 0.55

system.time(postA_1to500<-lapply(1:500, simulation_FUZZY,"B","0.55"))
system.time(postA_501to1000<-lapply(501:1000, simulation_FUZZY,"B","0.55"))

# Results
  Results1<-lapply(postB_1to500, performance_sample_FUZZY,"B","0.55")
  Results2<-lapply(postB_501to1000, performance_sample_FUZZY,"B","0.55")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results[,1:6], 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)
  sqrt(mean(Results[['c_abs_map']]^2)) # RMSE cutoff location (MAP estimate)
  sqrt(mean(Results[['j_abs_map']]^2)) # RMSE compliance rate (MAP estimate)


# LLR
LLR_res=LLR_performance_FUZZY(1000,'B','0.55')
apply(LLR_res, 2,mean)
apply(LLR_res[1:3], 2,median)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. 
sqrt(mean(LLR_res[['j_abs']]^2)) # RMSE compliance rate 


# jump = 0.3

system.time(postB_1to500<-lapply(1:500, simulation_FUZZY,"B","0.3"))
system.time(postB_501to1000<-lapply(501:1000, simulation_FUZZY,"B","0.3"))

# Results
  Results1<-lapply(postB_1to500, performance_sample_FUZZY,"B","0.3")
  Results2<-lapply(postB_501to1000, performance_sample_FUZZY,"B","0.3")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results[,1:6], 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)
  sqrt(mean(Results[['c_abs_map']]^2)) # RMSE cutoff location (MAP estimate)
  sqrt(mean(Results[['j_abs_map']]^2)) # RMSE compliance rate (MAP estimate)


# LLR
LLR_res=LLR_performance_FUZZY(1000,'B','0.3')
apply(LLR_res, 2,mean)
apply(LLR_res[1:3], 2,median)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. 
sqrt(mean(LLR_res[['j_abs']]^2)) # RMSE compliance rate 


# Function C #

# jump = 0.55

system.time(postC_1to500<-lapply(1001:1500, simulation_FUZZY,"C","0.55"))
system.time(postC_501to1000<-lapply(1501:2000, simulation_FUZZY,"C","0.55"))

# Results
  Results1<-lapply(postC_1to500, performance_sample_FUZZY,"C","0.55")
  Results2<-lapply(postC_501to1000, performance_sample_FUZZY,"C","0.55")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results[,1:6], 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)
  sqrt(mean(Results[['c_abs_map']]^2)) # RMSE cutoff location (MAP estimate)
  sqrt(mean(Results[['j_abs_map']]^2)) # RMSE compliance rate (MAP estimate)


# LLR
LLR_res=LLR_performance_FUZZY(1000,'C','0.55',1000)
apply(LLR_res, 2,mean)
apply(LLR_res[1:3], 2,median)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. 
sqrt(mean(LLR_res[['j_abs']]^2)) # RMSE compliance rate 


# jump = 0.3

system.time(postC_1to500<-lapply(1001:1500, simulation_FUZZY,"C","0.3"))
system.time(postC_501to1000<-lapply(1501:2000, simulation_FUZZY,"C","0.3"))

# Results
  Results1<-lapply(postC_1to500, performance_sample_FUZZY,"C","0.3")
  Results2<-lapply(postC_501to1000, performance_sample_FUZZY,"C","0.3")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results[,1:6], 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)
  sqrt(mean(Results[['c_abs_map']]^2)) # RMSE cutoff location (MAP estimate)
  sqrt(mean(Results[['j_abs_map']]^2)) # RMSE compliance rate (MAP estimate)


# LLR
LLR_res=LLR_performance_FUZZY(1000,'C','0.3',1000)
apply(LLR_res, 2,mean)
apply(LLR_res[1:3], 2,median)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. 
sqrt(mean(LLR_res[['j_abs']]^2)) # RMSE compliance rate 

## Sharp design ##
# A single simulation takes around 45 seconds. In the following code we split simulations
# into two batches of 500.


# Function A #

system.time(postA_1to500<-lapply(2001:2500, simulation_SHARP,"A"))
system.time(postA_501to1000<-lapply(2501:3000, simulation_SHARP,"A"))

  # Results
  Results1<-lapply(postA_1to500, performance_sample_SHARP,"A","LoTTA")
  Results2<-lapply(postA_501to1000, performance_sample_SHARP,"A","LoTTA")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)


# LLR
LLR_res=LLR_performance_SHARP(1000,'A',2000)
apply(LLR_res, 2,mean)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. (MAP estimate)


# Function B #

system.time(postB_1to500<-lapply(1:500, simulation_SHARP,"B"))
system.time(postB_501to1000<-lapply(501:1000, simulation_SHARP,"B"))

# Results
  Results1<-lapply(postB_1to500, performance_sample_SHARP,"B","LoTTA")
  Results2<-lapply(postB_501to1000, performance_sample_SHARP,"B","LoTTA")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)


# LLR
LLR_res=LLR_performance_SHARP(1000,'B')
apply(LLR_res, 2,mean)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. (MAP estimate)


# Function C #

system.time(postC_1to500<-lapply(1001:1500, simulation_SHARP,"C"))
system.time(postC_501to1000<-lapply(1501:2000, simulation_SHARP,"C"))

# Results
  Results1<-lapply(postC_1to500, performance_sample_SHARP,"C","LoTTA")
  Results2<-lapply(postC_501to1000, performance_sample_SHARP,"C","LoTTA")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)


# LLR
LLR_res=LLR_performance_SHARP(1000,'C',1000)
apply(LLR_res, 2,mean)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. (MAP estimate)


## Appendix B ##

# Function lee #

# LoTTA
system.time(postlee_1to500<-lapply(1:500, simulation_SHARP,"lee"))
system.time(postlee_501to1000<-lapply(501:1000, simulation_SHARP,"lee"))

# Results
  Results1<-lapply(postlee_1to500, performance_sample_SHARP,"lee","LoTTA")
  Results2<-lapply(postlee_501to1000, performance_sample_SHARP,"lee","LoTTA")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)


# Cubic polynomial
system.time(postlee_1to500<-lapply(1:500, simulation_cubic_SHARP,"lee"))
system.time(postlee_501to1000<-lapply(501:1000, simulation_cubic_SHARP,"lee"))

# Results
  Results1<-lapply(postlee_1to500, performance_sample_SHARP,"lee","3poly")
  Results2<-lapply(postlee_501to1000, performance_sample_SHARP,"lee","3poly")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)
  

# LLR
LLR_res=LLR_performance_SHARP(1000,'lee')
apply(LLR_res, 2,mean)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. (MAP estimate)


# Function ludwig #

# LoTTA
system.time(postludwig_1to500<-lapply(1:500, simulation_SHARP,"ludwig"))
system.time(postludwig_501to1000<-lapply(501:1000, simulation_SHARP,"ludwig"))
  # Results
  Results1<-lapply(postludwig_1to500, performance_sample_SHARP,"ludwig","LoTTA")
  Results2<-lapply(postludwig_501to1000, performance_sample_SHARP,"ludwig","LoTTA")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  sqrt(mean(Results[['abs_err_map']]^2))

# Cubic polynomial 
system.time(postludwig_1to500<-lapply(1:500, simulation_cubic_SHARP,"ludwig"))
system.time(postludwig_501to1000<-lapply(501:1000, simulation_cubic_SHARP,"ludwig"))
  # Results 
  Results1<-lapply(postludwig_1to500, performance_sample_SHARP,"ludwig","3poly")
  Results2<-lapply(postludwig_501to1000, performance_sample_SHARP,"ludwig","3poly")
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  sqrt(mean(Results[['abs_err_map']]^2)) # RMSE tr.eff. (MAP estimate)
  

# LLR
LLR_res=LLR_performance_SHARP(1000,'ludwig')
apply(LLR_res, 2,mean)
sqrt(mean(LLR_res[['abs_err']]^2)) # RMSE tr.eff. (MAP estimate)

## Plots ##

# main manuscript

# outcome functions
x=seq(-1,1,0.001)
X_plot=sort(2*rbeta(500,2,4)-1)
Y_plot=funA_sample(X_plot)
Fun=funA(x)
plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=Fun)

(plA=ggscatter(plot.df, x = "X", y = "Y",
               color = "gray",
               size = 3, alpha = 0.6)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1)+labs(title="Function A",x='',y='')+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.y=element_text(face="italic"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))

Y_plot=funB_sample(X_plot)
Fun=funB(x)
plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=Fun)
(plB=ggscatter(plot.df, x = "X", y = "Y",
               color = "gray",
               size = 3, alpha = 0.6)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1)+labs(title="Function B",x='',y='')+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.y=element_text(face="italic"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))

Y_plot=funC_sample(X_plot)
Fun=funC(x)
plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=Fun)
(plC=ggscatter(plot.df, x = "X", y = "Y",
               color = "gray",
               size = 3, alpha = 0.6)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1)+labs(title="Function C",x='',y='')+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.y=element_text(face="italic"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))

ggarrange(plA,plB,plC,ncol=3)

# treatment functions

Fun=fun_prob30(x)
plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=Fun)
(pl30=ggplot(lineF.df, x="X",y="Y")+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1)+ylim(0,0.8)+labs(title="Treatment probabilty function 2",x='',y='')+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.y=element_text(face="italic"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))

Fun=fun_prob55(x)
plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=Fun)
(pl50=ggplot(lineF.df, x="X",y="Y")+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1)+ylim(0,0.8)+labs(title="Treatment probabilty function 1",x='',y='')+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.y=element_text(face="italic"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))

ggarrange(pl50,pl30,ncol=2)

# appendix B1

x=seq(-1,1,0.001)
X_plot=sort(2*rbeta(500,2,4)-1)
Y_plot=lee_sample(X_plot)
Fun=lee(x)
plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=Fun)

(pllee=ggscatter(plot.df, x = "X", y = "Y",
               color = "gray",
               size = 3, alpha = 0.6)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1)+labs(title="Lee function",x='',y='')+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.y=element_text(face="italic"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))

Y_plot=ludwig_sample(X_plot)
Fun=ludwig(x)
plot.df=data.frame('X'=X_plot,'Y'=Y_plot)
lineF.df=data.frame(X=x,Y=Fun)
(plludwig=ggscatter(plot.df, x = "X", y = "Y",
               color = "gray",
               size = 3, alpha = 0.6)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_point(aes(x=X,y=Y),data=lineF.df,size=0.1)+labs(title="Ludwig function",x='',y='')+theme_classic(base_size = 14)+theme(text = element_text(family = "serif"),axis.title.y=element_text(face="italic"),axis.title.x=element_text(face="italic"),axis.text = element_text(family="sans"),plot.title = element_text(hjust = 0.5),legend.text =element_text(size = 14)))


ggarrange(pllee,plludwig,ncol=2)


#---------------------------#
#  Section 5 & Appendix B.2 #
#---------------------------#

# In the following simulations the outcome functions do not matter for the cutoff detection as only treatment models are used
# However, for each outcome function different seeds are used, therefore for the sake of precision we repeat cutoff detection
# with the corresponding seeds.

# In the LoTTA treatment-only model one simulations takes around 45 seconds. In the simple model with two constant functions one simulation takes around 7 seconds. 


# Function A 

# LoTTA treatement-only model

  # jump = 0.55 
  
  system.time(postA_1to500<-lapply(2001:2500, simulation_treatment,"0.55"))
  system.time(postA_501to1000<-lapply(2501:3000, simulation_treatment,"0.55"))
    
    # Results for the cutoff detection

  Results1<-lapply(postA_1to500, performance_cutoff_sample_FUZZY)
  Results2<-lapply(postA_501to1000, performance_cutoff_sample_FUZZY)
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"A","0.55")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))

  # jump = 0.30

  system.time(postA_1to500<-lapply(2001:1500, simulation_treatment,"0.3"))
  system.time(postA_501to1000<-lapply(2501:3000, simulation_treatment,"0.3"))
    
    # Results for the cutoff detection

  Results1<-lapply(postA_1to500, performance_cutoff_sample_FUZZY)
  Results2<-lapply(postA_501to1000, performance_cutoff_sample_FUZZY)
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"A","0.3")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))


# Two constant functions
  
  # jump = 0.55 
  
  system.time(postA_1to1000_s<-lapply(2001:3000, simulation_treatment_simple,"0.55"))
    
    # Results for the cutoff detection (not in the manuscript)

  Results<-lapply(postA_1to1000_s, performance_cutoff_sample_FUZZY)
  Results=do.call(rbind.data.frame, Results)
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"A","0.55")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))

  # jump = 0.30

  system.time(postA_1to1000_s<-lapply(2001:3000, simulation_treatment_simple,"0.3"))
    
    # Results for the cutoff detection (not in the manuscript)

  Results<-lapply(postA_1to1000_s, performance_cutoff_sample_FUZZY)
  Results=do.call(rbind.data.frame, Results)
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"A","0.3")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))



# Function B 

# LoTTA treatement-only model

  # jump = 0.55 
  
  system.time(postB_1to500<-lapply(1:500, simulation_treatment,"0.55"))
  system.time(postB_501to1000<-lapply(501:1000, simulation_treatment,"0.55"))
    
    # Results for the cutoff detection

  Results1<-lapply(postB_1to500, performance_cutoff_sample_FUZZY)
  Results2<-lapply(postB_501to1000, performance_cutoff_sample_FUZZY)
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"B","0.55")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))

  # jump = 0.30

  system.time(postB_1to500<-lapply(1:500, simulation_treatment,"0.3"))
  system.time(postB_501to1000<-lapply(501:1000, simulation_treatment,"0.3"))
    
    # Results for the cutoff detection

  Results1<-lapply(postB_1to500, performance_cutoff_sample_FUZZY)
  Results2<-lapply(postB_501to1000, performance_cutoff_sample_FUZZY)
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"B","0.3")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))


# Two constant functions
  
  # jump = 0.55 
  
  system.time(postB_1to1000_s<-lapply(1:1000, simulation_treatment_simple,"0.55"))
    
    # Results for the cutoff detection (not in the manuscript)

  Results<-lapply(postB_1to1000_s, performance_cutoff_sample_FUZZY)
  Results=do.call(rbind.data.frame, Results)
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"B","0.55")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))

  # jump = 0.30

  system.time(postB_1to1000_s<-lapply(1:1000, simulation_treatment_simple,"0.3"))
    
    # Results for the cutoff detection (not in the manuscript)

  Results<-lapply(postB_1to1000_s, performance_cutoff_sample_FUZZY)
  Results=do.call(rbind.data.frame, Results)
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"B","0.3")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))

# Function C 

# LoTTA treatement-only model

  # jump = 0.55 
  
  system.time(postC_1to500<-lapply(1001:1500, simulation_treatment,"0.55"))
  system.time(postC_501to1000<-lapply(1501:2000, simulation_treatment,"0.55"))
    
    # Results for the cutoff detection

  Results1<-lapply(postC_1to500, performance_cutoff_sample_FUZZY)
  Results2<-lapply(postC_501to1000, performance_cutoff_sample_FUZZY)
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"C","0.55")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))

  # jump = 0.30

  system.time(postC_1to500<-lapply(1001:1500, simulation_treatment,"0.3"))
  system.time(postC_501to1000<-lapply(1501:2000, simulation_treatment,"0.3"))
    
    # Results for the cutoff detection

  Results1<-lapply(postC_1to500, performance_cutoff_sample_FUZZY)
  Results2<-lapply(postC_501to1000, performance_cutoff_sample_FUZZY)
  Results1=do.call(rbind.data.frame, Results1)
  Results2=do.call(rbind.data.frame, Results2)
  Results=rbind(Results1,Results2) 
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"C","0.3")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))


# Two constant functions
  
  # jump = 0.55 
  
  system.time(postC_1to1000_s<-lapply(1001:2000, simulation_treatment_simple,"0.55"))
    
    # Results for the cutoff detection (not in the manuscript)

  Results<-lapply(postC_1to1000_s, performance_cutoff_sample_FUZZY)
  Results=do.call(rbind.data.frame, Results)
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"C","0.55")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))

  # jump = 0.30

  system.time(postC_1to1000_s<-lapply(2001:3000, simulation_treatment_simple,"0.3"))
    
    # Results for the cutoff detection (not in the manuscript)

  Results<-lapply(postC_1to1000_s, performance_cutoff_sample_FUZZY)
  Results=do.call(rbind.data.frame, Results)
  apply(Results, 2,mean)
  apply(Results, 2,median)
  sqrt(mean(Results[['abs_err_map']]^2)) 
    
    # Results for the treatment effect estimation with the plug-in estimator
  
  Results_plugin<-lapply(1:1000, plugin_estimator,Results,"C","0.3")
  Results_plugin=do.call(rbind.data.frame, Results_plugin)
  apply(Results_plugin, 2,mean)
  sqrt(mean(Results_plugin[['abs_err']]^2))


#----------------------------#
#  Appendix E.1              #
#----------------------------#





