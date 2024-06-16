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

data <- read.dta("/Users/julia/Downloads/CKT_2023_Chemo.dta", warn.missing.labels = FALSE)
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
ggscatter(plot.df, x = "X", y = "Y",
          color = "gray",
          size = 5, alpha = 0.6)+labs(title="Prob. of receiving chemotherapy",
                                      x ='Oncotype score', y = expression(paste("Proportion ","T=1")))+theme_classic(base_size = 13)+theme(text = element_text(family = "serif"),axis.text=element_text(family = "sans"))

# Scatter plot of the outcome data
Y_plot=Bin_data(Y,X,c,binsize)[[1]]
plot.df=data.frame(X=X_plot,Y=Y_plot)
ggscatter(plot.df, x = "X", y = "Y",
              color = "gray",
              size = 5, alpha = 0.6)+labs(title="Prob. of cancer recurrence",
                                          x ='Oncotype score', y = expression(paste("Proportion ","Y=1")))+theme_classic(base_size = 13)+theme(text = element_text(family = "serif"),axis.text=element_text(family = "sans"))

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

# Fit constant model
datchemo_T=list(N=length(X),x=X,t=T,jlb=0.2,prob=prob,start=start,nc=nc)
initcCHEMO=list(ct=10,al=0.5,j=0.3,.RNG.seed=1,.RNG.name="base::Mersenne-Twister")
param_c=c('ct')
datchemo_c<- run.jags('cutoff_initial_DIS.txt' ,data=datchemo_T,monitor=param_c,inits=list(initcCHEMO),burnin = 1000,sample=1000,adapt = 1000,n.chains = 1,method = 'parallel')
Ct_start=as.numeric(combine.mcmc(datchemo_c))
init1=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,100,nc)
init2=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,200,nc)
init3=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,300,nc)
init4=Initial_treatment_DIS(X,T,Ct_start,lb,ubr,ubl,start,prob,400,nc)
param_cj=c('c','j')
dat1=list(N=length(X),x=X,t=T,jlb=0.1,prob=prob,start=start,nc=nc,lb=lb,ubrt=ubr,ublt=ubl)
Chemo_treatment2<- run.jags('treatment_model_discrete.txt', data=dat1,monitor=param_cj,inits =list(init1,init2,init3,init4),burnin = 30000,sample=25000,adapt = 500, method='parallel',n.chains = 4)



# Fit treatment model with jlb=0.2