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
cloglog_IFR[k] ~ dnorm(theta[k], inv.var_tau);
cloglog_infectionrate[k] ~ dnorm(beta[k], inv.var_sig) ;
theta[k] <- alpha + omega * (ses[k] - 39.25029)  + ita1 * (age[k]-2)*(age[k]-3)*(age[k]-4) 
             + ita2 * (age[k]-1)*(age[k]-3)*(age[k]-4)
             + ita3 * (age[k]-1)*(age[k]-2)*(age[k]-4) + ita4 * (age[k]-2)*(age[k]-3)*(age[k]-4)* ses[k]
             + ita5 * (age[k]-1)*(age[k]-3)*(age[k]-4)* ses[k]
             + ita6 * (age[k]-1)*(age[k]-2)*(age[k]-4) * ses[k]
beta[k] <- balpha + bomega * (ses[k] - 39.25029)  + bita1 * (age[k]-2)*(age[k]-3)*(age[k]-4) 
             + bita2 * (age[k]-1)*(age[k]-3)*(age[k]-4)
             + bita3 * (age[k]-1)*(age[k]-2)*(age[k]-4) + bita4 * (age[k]-2)*(age[k]-3)*(age[k]-4)* ses[k]
             + bita5 * (age[k]-1)*(age[k]-3)*(age[k]-4)* ses[k]
             + bita6 * (age[k]-1)*(age[k]-2)*(age[k]-4) * ses[k]
IFR[k] <- icloglog(cloglog_IFR[k])
cloglog(infectionrate[k]) <- cloglog_infectionrate[k];

confirmed_cases[k] ~ dbin(infectionrate[k],tests[k]);
cases[k] ~ dbin(infectionrate[k], population[k]);
deaths[k] ~ dbin(IFR[k], cases[k]);
}

#define priors, here icloglog is the inverse link function g^-1
for (k in 1:K){ phi[k]<-1;}


#icloglog_alpha ~ dunif(0, 1); 
omega ~ dnorm(0, 1);
alpha ~ dnorm(0, 1);
ita1 ~ dnorm(0,1);
ita2 ~ dnorm(0,1);
ita3 ~ dnorm(0,1);
ita4 ~ dnorm(0,1);
ita5 ~ dnorm(0,1);
ita6 ~ dnorm(0,1);
bomega ~ dnorm(0, 1);
balpha ~ dnorm(0, 1);
bita1 ~ dnorm(0,1);
bita2 ~ dnorm(0,1);
bita3 ~ dnorm(0,1);
bita4 ~ dnorm(0,1);
bita5 ~ dnorm(0,1);
bita6 ~ dnorm(0,1);
#icloglog_beta0 ~ dunif(0, 1);

#icloglog_ita1 ~ dunif(0, 1);
#icloglog_ita2 ~ dunif(0, 1);
#icloglog_ita3 ~ dunif(0, 1);
#icloglog_ita4 ~ dunif(0, 1);
#icloglog_ita5 ~ dunif(0, 1);
#icloglog_ita6 ~ dunif(0, 1);


#alpha <- log(-log(1-icloglog_alpha));
#beta0 <- log(-log(1-icloglog_beta0));
#ita1 <- log(-log(1-icloglog_ita1));
#ita2 <- log(-log(1-icloglog_ita2));
#ita3 <- log(-log(1-icloglog_ita3));
#ita4 <- log(-log(1-icloglog_ita4));
#ita5 <- log(-log(1-icloglog_ita5));
#ita6 <- log(-log(1-icloglog_ita6));



inv.var_sig   <- (1/sd_sig)^2 ;
sd_sig     ~ dnorm(0, 1/1) T(0,);
inv.var_tau   <- (1/sd_tau)^2 ;
sd_tau     ~ dnorm(0, 1/0.01) T(0,);
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
  params <- c("phi", "alpha", "beta","omega",  "IFR", "infectionrate", "sd_sig",   "sd_tau","ita1"
              ,"ita2", "ita3","ita4","ita5","ita6","bomega","balpha", "bita1","bita2","bita3","bita4","bita5","bita6")		
  
  samps_ignore <- coda.samples(jags.m_ignore, params[-which.max(params=="gamma")], n.iter=MCMCiter,  n.burnin = MCMCiter*0.2, thin=50)
  median_ignore <- summary(samps_ignore)$quantiles[,c(3)]
  
  HDI<-hdi(samps_ignore)
  QQ_ignore<-t(rbind(HDI[1,],median_ignore,HDI[2,]))
  QQA<-QQ_ignore
  goodQ<-QQA[c( "beta",   "sd_tau", "sd_sig" ),]
  
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

MCMCtrace(results$samps_ignore, params= c("theta","beta0"), priors=cbind(cloglog(runif(MCMCiter_sim)),cloglog(runif(MCMCiter_sim))) , main_den = c(		   
  TeX("Density $\\theta$"),
  TeX("Density $\\beta$")),
  
  main_tr = c(					   					                         
    TeX("Trace $\\theta$"),
    TeX("Trace $\\beta$")),
  ,
  filename= "MCMC_theta_beta0_nonprefer.pdf")


MCMCtrace(results$samps_ignore, params= c("ita1","ita2","ita3"), priors=cbind(cloglog(runif(1)),cloglog(runif(1)),cloglog(runif(1))) , main_den = c(		   
  TeX("Density $\\ita1$"),
  TeX("Density $\\ita2$"),
  TeX("Density $\\ita2$")),
  
  main_tr = c(					   					                         
    TeX("Trace $\\ita1$"),
    TeX("Trace $\\ita2$"),
    TeX("Trace $\\ita2$"))
  ,
  filename= "MCMC_theta_ita_nonprefer.pdf")


IR <- cbind(results$QQA[(n+2):(2*n+1),],data$ses,data$Age_group)
colnames(IR)<-c('IRl','IRm','IRu','ses','age')

### Plot infection rate
IR %>% as_tibble() %>% ggplot(aes(x=ses,y=IRm)) +geom_point()+geom_errorbar(aes(ymin=IRl,ymax=IRu))

IR %>% as_tibble() %>% filter() %>% ggplot(aes(x=ses,y=IRm,color=factor(age))) +geom_point(position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=IRl,ymax=IRu))+facet_wrap(~age)+geom_smooth(method=lm)+
  stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~`,`~")))
