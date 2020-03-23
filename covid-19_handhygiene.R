# =============================================================================== #
# Hand hygiene covid-19 model 
# =============================================================================== #
rm(list=ls())
# Packages
library(ggplot2);library(tidyr);library(dplyr)

# Set your working directory here:
setwd("/home/tmp/Downloads/owncloud/PHD/Utrecht/JuliusCenter/Combacte/Coding/R/")
# Load functions
source("covid-19_handhygiene_functions.R")


###############################
# Sensitivity analysis for different secondary attack rates
epsMat <- sens.transPerContact.SAR(dc=seq(1,30)/60,infPeriod=12,c_mean=1,f_rate=10,
                                   sar_vec=c(0.15, 0.25, 0.5, 0.75))
plot.fun(epsMat, 
         x.title = "Half-life of probability of persistence (minutes)", 
         y.title = "Probability of transmission per contact",
         legend.title = "SAR", 
         legend.labels = c("15%","25%","50%","75%"))

# Sensitivity analysis for different rates of hand contamination events
# SAR = 0.5 fixed
epsMatC <- sens.transPerContact.cont(dc=seq(1,30)/60,infPeriod=12,c_mean=1,SAR=0.5,
                                     f_rate=10, c_vec=c(1/4,1/2,1,2,4,8,16))
plot.fun(epsMatC, 
         x.title = "Half-life of probability of persistence (minutes)", 
         y.title = "Probability of transmission per contact",
         legend.title = "Mean time between \nhand contamination \nevents", 
         legend.labels = c("15 min","30 min","1 hour","2 hours","4 hours","8 hours","16 hours"))


# =============================================================================== #
# how transmission prob in a 2 person household changes with hand hygiene 
# frequency as a function of duration of contamination (but where parameters are 
# constrained to give the same secondary attack rate in the absence of hand hygiene)
# =============================================================================== #
# Regular fixed hand washing
dataInf <- sar.hw(dc=seq(1,30)/60, 
                  hw=c(5/60,15/60,0.5,1,2,4,8,16), 
                  f_rate=10, 
                  infPeriod=12, 
                  c_mean=1, 
                  SAR.nohw=0.5,
                  HW.opt=1)

plot.fun(dataInf, 
         x.title = "Half-life of probability of persistence (minutes)", 
         y.title = "Secondary attack rate",
         legend.title = "Hand hygiene frequency \nEvery", 
         legend.labels = c("5 min", "15 min", "30 min", "1 hour", "2 hours", "4 hours", "8 hours", "16 hours"))


# Event-prompted hand washing
dataInf2 <- sar.hw(dc=seq(1,30)/60, 
                  hw=c(1/60,5/60,10/60,15/60,30/60), 
                  f_rate=10, 
                  infPeriod=12, 
                  c_mean=1, 
                  SAR.nohw=0.5,
                  HW.opt=1)

plot.fun(dataInf2, 
         x.title = "Half-life of probability of persistence (minutes)", 
         y.title = "Secondary attack rate",
         legend.title = "Delay between \nhand contamination \nand hand washing", 
         legend.labels = c("1 min","5 min","10 min", "15 min", "30 min"))




# =============================================================================== #
# plot of how many transmission events (in a population of 100 two person 
# households say) result from hands that have been contaminated for at least x 
# minutes in the absence of hand hygiene (so x axis is minutes, y axis is 
# proportion of transmission events)
# =============================================================================== #
epsilon <- 0.001 # Probability of transmission per contact
num <- 100 # Number of households
dc <- c(1,3,5,10,15)/60 # Half life of probability of persistence
propTransMat <- NULL
maxExp <- 6
xaxis <- seq(3/60,maxExp,by=3/60) # Discretization of time
for(d in dc){
  pinf <- rep(0,num)
  texp <- v <- last_tc <- list()
  for(h in 1:num){
    sim <- household.fun(t_end=24*infPeriod, HW=0, f_rate=11, eps=epsilon, halflife=d)
    pinf[h] <- sim$p_inf
    texp[[h]] <- sim$t_exp
    v[[h]] <- sim$v
    last_tc[[h]] <- sim$last_tc
  }
  v <- lapply(v, function(x) x[which(x!=0)])
  
  # maxExp <- max(c(maxExp, unlist(lapply(texp, function(x) max(round(x,2))))))
  propTrans <- rep(0, length(xaxis))
  for(i in 1:length(xaxis)){
    temp <- NULL
    for(h in 1:num){
      temp<- c(temp, 1-exp(-sum(epsilon*v[[h]][which(texp[[h]]<xaxis[i])])))
    }
    propTrans[i] <- sum(temp)/num
  }
  propTransMat <- cbind(propTransMat, propTrans)
}

propTransMat <- as.data.frame(propTransMat)
colnames(propTransMat) <- 1:ncol(propTransMat)
rownames(propTransMat) <- 1:nrow(propTransMat)
propTransMat %>% gather() %>% group_by(key) %>% 
  mutate(x=1:n()) %>%
  ggplot(aes(x=x, y=value,group=key,color=key)) + 
  geom_line() +
  xlab("Time of contaminated hands after contamination event (minutes)") + 
  ylab("Proportion of transmission events") + 
  # scale_x_discrete(breaks=xaxis, labels=as.character(xaxis)) + 
  scale_color_discrete(name="Half life of probability \nof persistence (minutes)",
                       labels=as.character(dc*60))

# =============================================================================== #
# exploration of whether distribution of duration of carriage duration matters. 
# i.e. are results the same if duration is constant or exponential 
# (with the same mean) - or heavy-tailed 
# =============================================================================== #
dc <- seq(1,30)/60
p <- q <- r <- rep(0, length(dc))
for(i in 1:length(dc)){
  p[i] <- mean(sapply(1:num, function(x)household.fun(t_end=24*infPeriod, HW=0, CT=1, c_mean=1,f_rate=11, eps=0.001, halflife=dc[i])$p_inf))
  q[i] <- mean(sapply(1:num, function(x)household.fun(t_end=24*infPeriod, HW=0, CT=2, c_mean=1,f_rate=11, eps=0.001, halflife=dc[i])$p_inf))
  r[i] <- mean(sapply(1:num, function(x)household.fun(t_end=24*infPeriod, HW=0, CT=3, c_mean=1,c_sd=1, f_rate=11, eps=0.001, halflife=dc[i])$p_inf))
}

dataSens <- as.data.frame(cbind(p,q,r)) 
colnames(dataSens) <- 1:ncol(dataSens)
dataSens %>% gather() %>% group_by(key) %>% 
  mutate(x=1:n()) %>%
  ggplot(aes(x=x, y=value,group=key,color=key)) + 
  geom_line() +
  xlab("Half life of probability \n of persistence (minutes)") + 
  ylab("Proportion of transmission events") + 
  scale_color_discrete(name="Distribution",
                       labels=c("Uniform","Exponential","Lognormal"))

