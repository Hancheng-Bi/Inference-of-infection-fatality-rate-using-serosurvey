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

### read in data from desktop
data <- read_csv("data_fig5C.csv")
View(data)
data<- data[1:4,]
n = length(data$Positive)

### JAGS model###
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

# Define priors
for (k in 1:K){ phi[k]<-1;}
theta ~ dnorm(0, 1);
beta ~ dnorm(0, 1);
sig     ~ dnorm(0, 1/1) T(0,);
tau     ~ dnorm(0, 1/1) T(0,);
}'
#    
cat(model, file="JAGS_nonpreferC.txt")
############
### END of JAGS model for "M_{0}" ###
############

illustrative <- function(observed_data, MCMCiter=500000){
  
  ###
  datalist <- list(K=length(observed_data$Positive), confirmed_cases=unlist(observed_data$Positive), deaths=observed_data$death, population=observed_data$Pop, tests=observed_data$Sample_size)
  
  jags.m_ignore <- jags.model(file = "JAGS_nonpreferC.txt", data = datalist, n.chains = 5, n.adapt = 50000, inits= NULL)
  
  ######
  params <- c("phi", "theta", "beta",  "IFR", "infectionrate", "sig",   "tau")		
  
  samps_ignore <- coda.samples(jags.m_ignore, params[-which.max(params=="gamma")], n.iter=MCMCiter,  n.burnin = MCMCiter*0.2, thin=50)
  median_ignore <- summary(samps_ignore)$quantiles[,c(3)]
  
  HDI<-hdi(samps_ignore)
  QQ_ignore<-t(rbind(HDI[1,],median_ignore,HDI[2,]))
  QQA<-QQ_ignore
  goodQ<-QQA[c("theta", "beta",   "tau", "sig" ),]
  
  return(list(goodQ = goodQ, QQA=QQA, samps_ignore = samps_ignore))
}
############
### end of illustrative function ###
############

### Run MCMC
## "Each model is fit using JAGS (just another Gibbs' sampler) \citep{kruschke2014doing}, with 5 independent chains, each with 500,000 draws (20\% burnin, thinning of 50)." ###  

MCMCiter_sim <- 500000
results  <- illustrative(data,  MCMCiter = MCMCiter_sim)
n_sample <- MCMCiter_sim/50
IFR_chain <- rbind(results$samps_ignore[1][1:n_sample,1:n][[1]],results$samps_ignore[2][1:n_sample,1:n][[1]],
                   results$samps_ignore[3][1:n_sample,1:n][[1]],results$samps_ignore[4][1:n_sample,1:n][[1]],
                   results$samps_ignore[5][1:n_sample,1:n][[1]])

colnames(IFR_chain) <- c('IFR1','IFR2','IFR3','IFR4')




### read in data from desktop
data2 <- read_csv("data_fig5C.csv")
View(data2)
data2<- data2[5:8,]
n = length(data$Positive)

### JAGS model for "M_{0}" ###
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

for (k in 1:K){ phi[k]<-1;}

theta ~ dnorm(0, 1);
beta ~ dnorm(0, 1);
sig     ~ dnorm(0, 1/1) T(0,);
tau     ~ dnorm(0, 1/1) T(0,);
}'
#    
cat(model, file="JAGS_nonpreferC.txt")
############
### END of JAGS model for "M_{0}" ###
############

illustrative <- function(observed_data, MCMCiter=500000){
  
  ###
  datalist <- list(K=length(observed_data$Positive), confirmed_cases=unlist(observed_data$Positive), deaths=observed_data$death, population=observed_data$Pop, tests=observed_data$Sample_size)
  # ignore is when model assumes gamma = 0
  jags.m_ignore <- jags.model(file = "JAGS_nonpreferC.txt", data = datalist, n.chains = 5, n.adapt = 50000, inits= NULL)
  
  ######
  params <- c("phi", "theta", "beta",  "IFR", "infectionrate", "sig",   "tau")		
  
  samps_ignore <- coda.samples(jags.m_ignore, params[-which.max(params=="gamma")], n.iter=MCMCiter,  n.burnin = MCMCiter*0.2, thin=5)
  median_ignore <- summary(samps_ignore)$quantiles[,c(3)]
  
  HDI<-hdi(samps_ignore)
  QQA<-t(rbind(HDI[1,],median_ignore,HDI[2,]))
  
  goodQ<-QQA[c("theta", "beta",   "tau", "sig" ),]
  
  return(list(goodQ = goodQ, QQA=QQA, samps_ignore = samps_ignore))
}
############
### end of illustrative function ###
############

### Run MCMC
## "Each model is fit using JAGS (just another Gibbs' sampler) 
## \citep{kruschke2014doing}, with 5 independent chains, each with 500,000 draws (20\% burnin, thinning of 50)." ###  

MCMCiter_sim <- 500000
results  <- illustrative(data2,  MCMCiter = MCMCiter_sim)
n_sample <- MCMCiter_sim/50
IFR_chain1 <- rbind(results$samps_ignore[1][1:n_sample,1:n][[1]],results$samps_ignore[2][1:n_sample,1:n][[1]],
                    results$samps_ignore[3][1:n_sample,1:n][[1]],results$samps_ignore[4][1:n_sample,1:n][[1]],
                    results$samps_ignore[5][1:n_sample,1:n][[1]])

colnames(IFR_chain1) <- c('IFR5','IFR6','IFR7','IFR8')

IFR_chain2<-cbind(IFR_chain,IFR_chain1)


## Compute ratio of posterior of high ses over low ses
ratio1<-IFR_chain2[,1]/IFR_chain2[,5]
ratio2 <- IFR_chain2[,2]/IFR_chain2[,6]
ratio3<-IFR_chain2[,3]/IFR_chain2[,7]
ratio4 <- IFR_chain2[,4]/IFR_chain2[,8]
age <- c('1','2','3','4')
ratem <- c(median(ratio1),median(ratio2),median(ratio3),median(ratio4))
ratel <- c(quantile(ratio1,0.05/2),quantile(ratio2,0.05/2),quantile(ratio3,0.05/2),quantile(ratio4,0.05/2))
rateu <- c(quantile(ratio1,1-0.05/2),quantile(ratio2,1-0.05/2),quantile(ratio3,1-0.05/2),quantile(ratio4,1-0.05/2))
IFR_rate <- cbind(age,ratem,ratel,rateu)
colnames(IFR_rate) <- c('age','ratem','ratel','rateu')
write.csv(IFR_rate,file = 'IFR_rate.csv')
