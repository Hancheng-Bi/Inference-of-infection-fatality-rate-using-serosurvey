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
ilogit <- function(x) {exp(x)/(1+exp(x))}
### read in data from desktop
data <- read_csv("data_combined.csv")

View(data)
n = length(data$Positive)
ses.mean = mean(data_combined$ses)


### JAGS model ###
############
model<- '
model{ 
#likelihood
for (k in 1:K){
cloglog_IFR[k] ~ dnorm(theta[k], (tau)^2);
cloglog_infectionrate[k] ~ dnorm(beta[k], (sig)^2) ;
theta[k] <- alpha0 + alpha1 * ses[k]  + ita[age[k]] + gamma[age[k]] + epsilon

IFR[k] <- icloglog(cloglog_IFR[k])
cloglog(infectionrate[k]) <- cloglog_infectionrate[k];

confirmed_cases[k] ~ dbin(infectionrate[k],tests[k]);
cases[k] ~ dbin(infectionrate[k], population[k]);
deaths[k] ~ dbin(IFR[k], cases[k]);
}

#define priors, here icloglog is the inverse link function g^-1
for (k in 1:K){ phi[k]<-1;}



alpha0 ~ dnorm(0, 1);
alpha1 ~ dnorm(0, 1);
for (j in 1:4){
ita[j] ~ dnorm(0,1)；
gamma[j] ~ dnorm(0,1)
};
epsilon ~ dnorm(0,1)




sig     ~ dnorm(0, 1/1) T(0,);
tau     ~ dnorm(0, 1/1) T(0,);
}'
#    
cat(model, file="JAGS.txt")
############
### END of JAGS model  ###
############

illustrative <- function(observed_data, MCMCiter=5000){
  
  ###
  datalist <- list(K=length(observed_data$Positive), confirmed_cases=unlist(observed_data$Positive), deaths=observed_data$death, population=observed_data$Pop, tests=observed_data$Sample_size
                   , ses = observed_data$ses, age = observed_data$Age_group )
  # ignore is when model assumes gamma = 0
  jags.m_ignore <- jags.model(file = "JAGS.txt", data = datalist, n.chains = 5, n.adapt = 5000, inits= NULL)
  
  ######
  params <- c("phi", "alpha0", "beta",  "IFR", "infectionrate", "sig",   "tau","ita"
              "alpha1","gamma", "epsilon")  
  
  samps_ignore <- coda.samples(jags.m_ignore, params[-which.max(params=="gamma")], n.iter=MCMCiter,  n.burnin = MCMCiter*0.2, thin=50)
  median_ignore <- summary(samps_ignore)$quantiles[,c(3)]
  
  HDI<-hdi(samps_ignore)
  QQ_ignore<-t(rbind(HDI[1,],median_ignore,HDI[2,]))
  QQA<-QQ_ignore
  goodQ<-QQA[c( "beta",   "tau", "sig" ),]
  
  return(list(goodQ = goodQ, QQA=QQA, samps_ignore = samps_ignore))
}
############
### end of illustrative function ###
############

### Run MCMC
## "Each model is fit using JAGS (just another Gibbs' sampler) \citep{kruschke2014doing}, with 5 independent chains, each with 500,000 draws (20\% burnin, thinning of 50)." ###  

MCMCiter_sim <- 5000
results  <- illustrative(data,  MCMCiter = MCMCiter_sim)


## Find quatiles of IFR
IFR  <- results$QQA[1:n,]
IFRm <- results$QQA[1:n,2]
IFRl <- results$QQA[1:n,1]
IFRu <- results$QQA[1:n,3]
IFR <- cbind(data$Age_group,data$ses,IFR)
colnames(IFR)<-c('Age','Ses','IFRl','IFRm','IFRu')

write.csv(IFR,'IFR_nonprefer.csv')
