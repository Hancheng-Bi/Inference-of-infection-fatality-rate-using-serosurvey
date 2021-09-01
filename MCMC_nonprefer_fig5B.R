library(readr)
library(rjags)
library(xtable)
library(MCMCvis)
library(coda)
library(MCMCpack)
library(BiasedUrn)
library(rvest)
library(runjags)
library(MASS)
library(HDInterval)
library(gridGraphics)
library(gridExtra)
library(forestplot)
library(latex2exp)
library(sf)
library(lemon)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(ggspatial)
library(lubridate)


### functions:
invlog <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
cloglog <- function(x){log(-log(1-x))}
icloglog <- function(x){1 - exp(-exp(x))}

### read in data 
data <- read_csv("data_combined.csv")
View(data)
n = length(data$Positive)

### JAGS text for non-preferential model###
############
model<- '
model{	
#likelihood
for (k in 1:K){
cloglog_IFR[k] ~ dnorm(theta, (tau)^2);
cloglog_infectionrate[k] ~ dnorm(beta, (sig)^2) ;
cloglog(IFR[k]) <- cloglog_IFR[k];
cloglog(infectionrate[k]) <- cloglog_infectionrate[k];
confirmed_cases[k] ~ dbin(infectionrate[k],tests[k]);
cases[k] ~ dbin(infectionrate[k], population[k]);
deaths[k] ~ dbin(IFR[k], cases[k]);
}

#define priors
for (k in 1:K){ phi[k]<-1;}
theta ~ dnorm(0, 1);
beta ~ dnorm(0, 1);

sig     ~ dnorm(0, 1/1) T(0,);
tau     ~ dnorm(0, 1/1) T(0,);
}'
#    
cat(model, file="JAGS_nonprefer.txt")
############
### END of  model 
############

# Main function
illustrative <- function(observed_data, MCMCiter=50000){
  
  ###
  datalist <- list(K=length(observed_data$Positive), confirmed_cases=unlist(observed_data$Positive), deaths=observed_data$death, population=observed_data$Pop, tests=observed_data$Sample_size)
  # ignore is when model assumes gamma = 0
  jags.m_ignore <- jags.model(file = "JAGS_nonprefer.txt", data = datalist, n.chains = 5, n.adapt = 5000, inits= NULL)
  
  ######
  params <- c("phi", "theta", "beta",  "IFR", "infectionrate", "sig",   "tau")		
  
  samps_ignore <- coda.samples(jags.m_ignore, params[-which.max(params=="gamma")], n.iter=MCMCiter,  n.burnin = MCMCiter*0.2, thin=50)
  median_ignore <- summary(samps_ignore)$quantiles[,c(3)]
  
  HDI<-hdi(samps_ignore)
  QQ_A<-t(rbind(HDI[1,],median_ignore,HDI[2,]))
  goodQ<-QQA[c("theta", "beta",   "tau", "sig" ),]
  
  return(list(goodQ = goodQ, QQA=QQA, samps_ignore = samps_ignore))
}
############
### end of illustrative function ###
############

### Run MCMC
## "Each model is fit using JAGS , with 5 independent chains, each with 500,000 draws 
## (20\% burnin, thinning of 50)." ###  
MCMCiter_sim <- 50000
results5B  <- illustrative(data,  MCMCiter = MCMCiter_sim)
## Find quatiles of IFR
IFR  <- results5B$QQA[1:n,]
IFRm <- results5B$QQA[1:n,2]
IFRl <- results5B$QQA[1:n,1]
IFRu <- results5B$QQA[1:n,3]
IFR <- cbind(data$Age_group,data$ses,IFR)
colnames(IFR)<-c('Age','Ses','IFRl','IFRm','IFRu')
write.csv(IFR,'IFR_nonprefer.csv')

### Trace plot
MCMCtrace(results5B$samps_ignore, params= c("theta","beta"), priors=cbind(cloglog(runif(MCMCiter_sim)),cloglog(runif(MCMCiter_sim))) , main_den = c(		   
  TeX("Density $\\theta$"),
  TeX("Density $\\beta$")),
  
  main_tr = c(					   					                         
    TeX("Trace $\\theta$"),
    TeX("Trace $\\beta$")),
    ,
  filename= "MCMC_theta_beta_nonprefer.pdf")

## Find parameters of MCMC
parameter <-cbind(c('theta','beta','tau','sig'), results5B$goodQ)
colnames(parameter) <- c('para','lower','median','upper')
write.csv(parameter,'parameter.csv')
  
  
  
  
  
