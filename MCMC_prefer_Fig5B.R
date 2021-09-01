## This file plot Fig5B that includes preferential testing parameters
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
data <- read_csv("Desktop/clean_data.csv")
View(data)
n = length(data$Positive)

### JAGS model for preferential testing ###
############
model<- '
model{	
#likelihood
for (k in 1:K){
cloglog_IFR[k] ~ dnorm(theta, (tau)^2);
cloglog_infectionrate[k] ~ dnorm(beta, (sig)^2) ;
cloglog(IFR[k]) <- cloglog_IFR[k];
cloglog(infectionrate[k]) <- cloglog_infectionrate[k];

## confirmed_cases[k] ~ dhyper(cases[k], population[k]-cases[k], tests[k], phi[k]);
confirmed_cases[k] ~ dbin(1-(1-infectionrate[k])^phi[k],tests[k]);
cases[k] ~ dbin(infectionrate[k], population[k]);
deaths[k] ~ dbin(IFR[k], cases[k]);
}

#priors
for (k in 1:K){ phi[k] ~ dunif(1, 1+gamma);}
gamma  ~ dexp(0.5);

theta ~ dnorm(0, 1);
beta ~ dnorm(0, 1);
sig     ~ dnorm(0, 1/1) T(0,);
tau     ~ dnorm(0, 1/1) T(0,);
}'
#    
cat(model, file="JAGS_prefer.txt")
############
### JAGS model for preferential testing ###
############

illustrativep <- function(observed_data, MCMCiter=50000){
  
  ###
  datalist <- list(K=length(observed_data$Positive), confirmed_cases=unlist(observed_data$Positive), deaths=observed_data$death, population=observed_data$Pop, tests=observed_data$Sample_size)
  jags.m_unif <- jags.model(file = "JAGS_prefer.txt", data = datalist, n.chains = 5, n.adapt = 50000, inits= NULL)
  
  ######
  params <- c("phi", "theta", "beta",  "IFR", "infectionrate", "sig",  "gamma", "tau")		
  
  samps_unif <- coda.samples(jags.m_unif, params, n.iter=MCMCiter,  n.burnin = MCMCiter*0.2, thin=50)
  
  median_unif <- summary(samps_unif)$quantiles[,c(3)]
  
  HDI<-hdi(samps_unif)
  QQ_A<-t(rbind(HDI[1,],median_unif,HDI[2,]))
  goodQ<-QQA[c("theta", "beta",   "tau", "sig" , "gamma"),]
  
  return(list(goodQ = goodQ, QQA=QQA, samps_unif = samps_unif))
}
############
### end of illustrative function ###
############

### Run MCMCcitep{kruschke2014doing}, with 5 independent chains, 
## each with 500,000 draws (20\% burnin, thinning of 50)." ###  

MCMCiter_sim <- 50000
result  <- illustrativep(data,  MCMCiter = MCMCiter_sim)


## Find quatiles of IFR
IFR  <- result$QQA[1:n,]
IFRm <- result$QQA[1:n,2]
IFRl <- result$QQA[1:n,1]
IFRu <- result$QQA[1:n,3]
IFR <- cbind(data$Age_group,data$ses,IFR)
colnames(IFR) <- c('Age','Ses','IFRl','IFRm','IFRu')
write.csv(IFR,'IFR_prefered.csv')

