# =============================================================================== #
# Hand hygiene covid-19 model 
# =============================================================================== #
rm(list=ls())
# Packages
library(ggplot2);library(tidyr);library(dplyr); library(gridExtra)
library(grid);

# Set your working directory here:
dataPath <- "/home/thi.mui.pham/covid-19/handhygiene/modelling/code/"
plotPath <- "/home/thi.mui.pham/covid-19/handhygiene/modelling/figures/"

# Load functions
source(paste0(dataPath,"covid-19_handhygiene_functions.R"))

infPeriod <- 12
c_mean <- 1/4
f_rate <- 10
hl1 <- 5.4/60; hl2 <- 36.1/60
sarNoHW <- 0.1
seed <- 100

eps <- compute.eps(t_end=infPeriod,c_mean=c_mean,f_rate=f_rate,halflife=hl1,SAR=sarNoHW, seed=seed)
eps
test <- household.fun(t_end=infPeriod, 
                      CT=2, 
                      c_mean=c_mean, 
                      halflife=hl1, 
                      hw_mean = 7.5/60, 
                      eps=eps, seed=100)
test$p_inf

eps2 <- compute.eps(t_end=infPeriod,c_mean=c_mean,f_rate=f_rate,halflife=hl2,SAR=sarNoHW, seed=seed)
test2 <- household.fun(t_end=infPeriod, 
                       CT=2, 
                       c_mean=c_mean, 
                       halflife=hl2, 
                       hw_mean = 15/60, 
                       eps = eps2, seed=100)
test2$p_inf

summary(test$t_exp)
summary(test2$t_exp)

par(mfrow=c(1,2))
boxplot(test$t_exp)
boxplot(test2$t_exp)

sum(test$t_exp)
sum(test2$t_exp)

View(cbind(test$v, test2$v))

eps <- compute.eps(t_end=infPeriod,c_mean=c_mean,f_rate=f_rate,halflife=hl1,SAR=sarNoHW, seed=seed)
test <- household.fun(t_end=infPeriod, 
                      CT=2, 
                      c_mean=c_mean, 
                      halflife=hl1, 
                      HW=0,
                      eps=eps, seed=100)

eps2 <- compute.eps(t_end=infPeriod,c_mean=c_mean,f_rate=f_rate,halflife=hl2,SAR=sarNoHW, seed=seed)
test2 <- household.fun(t_end=infPeriod, 
                       CT=2, 
                       c_mean=c_mean, 
                       halflife=hl2, 
                       HW=0, 
                       eps = eps2, seed=100)
test$p_inf
test2$p_inf
sum(test$t_exp)
sum(test$t_exp)

hl3 <- 1/60
eps3 <- compute.eps(t_end=infPeriod,c_mean=c_mean,f_rate=f_rate,halflife=hl3,SAR=sarNoHW, seed=seed)
test3 <- household.fun(t_end=infPeriod, 
                       CT=2, 
                       c_mean=c_mean, 
                       halflife=hl3, 
                       hw_mean = 2/60, 
                       eps = eps3, seed=100)
test3$p_inf

summary(test$t_exp)
summary(test2$t_exp)
summary(test3$t_exp)

# =============================================================================== #
# What does model predict about relationship between hand contamination rate,
# half life of contamination, and relative risk (i.e. effect estimated from 
# meta-regression which is RR of infection per day associated with one hand 
# hygiene event)
# =============================================================================== #
infPeriod <- 0.5
sarNoHW <- 0.01
hl1 <- 5.4; hl2 <- 36.1
# Regular fixed hand washing
data <- df_red <- NULL
for(c in c(10/60,30/60,1,1.5,2,2.5,3,3.5,4,4.5,5)){
  dataInf <- sar.hw(dc=seq(1,60)/60,
                    hw=2, 
                    # hw=c(5/60,15/60,0.5,1,2,4,6), 
                    f_rate=10, 
                    infPeriod=infPeriod, 
                    c_mean=c, 
                    SAR.nohw=sarNoHW,
                    HW.opt=1, seed=1,it=10)
  df <- dataInf$dataInf
  data <- rbind(data, t(df))
  df_red <- rbind(df_red, t((sarNoHW-df)/sarNoHW))
}

df_red
x <- c(10/60,30/60,1,1.5,2,2.5,3,3.5,4,4.5,5)
y <- seq(1,60)/60

plot_ly(x=x,y=y,z=as.matrix(df_red),type="contour")

1-(1-0.04)^6

RRmean <- 9.588737e-01
RRlow <- 7.864217e-01 
RRup <- 1.122748

1-RRmean^6
1-RRlow^6
1-RRup^6


# Vary baseline probability of infeciton (sarNoHW)
data <- red <- NULL
for(s in c(0.01,seq(0.1,0.9,by=0.1))){
  dataInf <- sar.hw(dc=seq(1,60)/60,
                    hw=, 
                    # hw=c(5/60,15/60,0.5,1,2,4,6), 
                    f_rate=10, 
                    infPeriod=infPeriod, 
                    c_mean=4, 
                    SAR.nohw=s,
                    HW.opt=1, seed=100,it=10)
  df <- dataInf$dataInf
  data <- rbind(data, t(df))
  red <- rbind(red, t((s-df)/s))
}


red[,36]
red[,60]










#################
sim <- household.fun(t_end=24*1,HW=3,CT=2,t_delay=5,f_rate=10,eps=0.001,halflife=5/60,seed=12345)
length(sim$t_HW)
sim2 <- household.fun(t_end=24*12,HW=1,hw_mean=30/60,f_rate=10,eps=0.001,halflife=5/60,seed=12345)
length(sim2$t_HW)

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
    sim <- household.fun(t_end=infPeriod, HW=0, f_rate=11, eps=epsilon, halflife=d)
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


# ====================================================================== #
# Sensitivity analysis for different secondary attack rates
epsMat <- sens.transPerContact.SAR(dc=seq(1,30)/60,infPeriod=12,c_mean=1,f_rate=10,
                                   sar_vec=c(0.15, 0.25, 0.5, 0.75))
plot.fun(epsMat, 
         x.title = "Half-life of probability of persistence (minutes)", 
         y.title = "Probability of transmission per contact",
         legend.title = "SAR", 
         legend.labels = c("15%","25%","50%","75%"))

# ====================================================================== #
# Sensitivity analysis for different rates of hand contamination events
# SAR = 0.5 fixed
epsMatC <- sens.transPerContact.cont(dc=seq(1,30)/60,infPeriod=12,c_mean=1,SAR=0.5,
                                     f_rate=10, c_vec=c(1/4,1/2,1,2))
plot.fun(epsMatC, 
         x.title = "Half-life of probability of persistence (minutes)", 
         y.title = "Probability of transmission per contact",
         legend.title = "Mean time between \nhand contamination", 
         legend.labels = c("15 min","30 min","1 hour","2 hours"))

