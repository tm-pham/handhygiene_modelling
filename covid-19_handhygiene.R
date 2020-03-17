# =============================================================================== #
# Hand hygiene covid-19 model 
# =============================================================================== #
rm(list=ls())
# Packages
library(ggplot2)
library(tidyr)
library(dplyr)

# Generate times according to exponential distribution
# Exponentially distributed time-to-event until time t_end
exp.times <- function(t_end, rate){
  t <- t_temp <- NULL
  t0 <- t_max <- 0
  while(t_max< t_end){
    temp <- rexp(1, rate=rate)
    t_temp <- c(t_temp, temp)
    t_max <- t_max + temp
    t <-c(t, t_max) 
  }
  if(t[length(t)]>t_end) {
    t <- t[-length(t)]
    t_temp <- t_temp[-length(t_temp)]
  }
  return(list(t=t,t_temp=t_temp))
}

lognormal.times <- function(t_end, c_mean, c_sd){
  t <- t_temp <- NULL
  t0 <- t_max <- 0
  while(t_max< t_end){
    temp <- rlnorm(1, meanlog=c_mean, sdlog=c_sd)
    t_temp <- c(t_temp, temp)
    t_max <- t_max + temp
    t <-c(t, t_max) 
  }
  if(t[length(t)]>t_end) {
    t <- t[-length(t)]
    t_temp <- t_temp[-length(t_temp)]
  }
  return(list(t=t,t_temp=t_temp))
}

# =============================================================================== #
# Household function
# =============================================================================== #
# unit of time = hour
# c_rate = Rate at which hands of person is contaminated
# f_rate = Rate at which hands touch face (per unit of time)
# halflife = Time until probability of persistence is halved (per unit of time)
# f_times = Mean number of times at which hands touch face (per unit of time)
# t_end = Time period (default: 24 hours)
# eps = Probability of transmission per contact
# HW = Option for hand washing frequency (1=fixed regular times, 2=exponential)
# For HW=1: t_reg = regular, fixed times between hand washing
# For HW=2: hw_mean = mean for exponential distribution of hand washing times
# CT = Option for distribution of carriage duration
# =============================================================================== #

household.fun <- function(t_end=24*5, 
                          CT = 2, c_mean=1, c_sd=1, halflife=1, 
                          f_rate=2, eps=0.01,
                          HW = 1, hw_mean=10/60,
                          seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  # rate for exponential decay computed from halflife
  v_rate <- -log(0.5)/halflife
  
  # Contamination times
  if(CT==1){ # Fixed times (constant)
    t_c<- seq(0, t_end, by = c_mean)
  }else if(CT==2){ # Exponential distrbution
    t_c <- exp.times(t_end,1/c_mean)$t
  }else if(CT==3){ # Log-normal (heavy-tail)
    t_c <- lognormal.times(t_end,c_mean,c_sd)$t
  }
  
  
  # Hand washing times/events
  t_HW <-NULL
  if(HW==1){# Hand washing at fixed times
    t_HW <- seq(0, t_end, by = hw_mean)
  }else if(HW==2){# Hand washing at random times (Poisson)
    t_HW <- exp.times(t_end, 1/hw_mean)$t
  }else if(HW==3){# Event-prompted hand washing
    # TBD
  }
  
  # Face-touching events
  # For now: poisson distribution for events/exponential distribution for times
  t_f <- exp.times(t_end, f_rate)$t
  
  # Compute probability of persistence at face-touching events
  # For now: Exponential distribution 
  v <- rep(0, length(t_f)) # prob. of contamination at face touching events
  # total time of exposure, i.e. cumulative time between contamination event and face-touching 
  # with no hand washing event: 
  t_exp <- last_tc <- NULL
  for(i in 1:length(t_f)){
    ind <- which(t_c<=t_f[i])
    if(length(ind)>0){
      # Last contamination event/time
      last <- t_c[ind[length(ind)]] 
      if(HW==0){
        v[i] <- exp(-v_rate*(t_f[i]-last)) # No handwashing events
        t_exp <- c(t_exp, t_f[i]-last)
        last_tc <- c(last_tc,last)
      }
      else{
        # Was there a hand washing in between contamination and face-touching?
        ind_HW <- intersect(which(t_HW>=last), which(t_HW<t_f[i]))
        if(length(ind_HW)>0) v[i] <- 0 # yes
        else{
          v[i] <- exp(-v_rate*(t_f[i]-last)) # no
          t_exp <- c(t_exp, t_f[i]-last)
          last_tc <- c(last_tc,last)
        }
      }
    }
  }
  
  # Force of infection
  foi <- eps*v
  # Probability of infection
  p_inf <- 1-exp(-sum(foi))
  
  return(list(p_inf=p_inf, v=v, t_c=t_c, t_f=t_f, t_HW=t_HW, t_exp=t_exp, last_tc=last_tc))
}

# Calculate epsilon (prob. of transmission per contact) in scenario with no HW
# Fixed SAR
compute.eps <- function(t_end,c_mean,f_rate,halflife, SAR=0.5, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  t_c <- exp.times(t_end,1/c_mean)$t
  t_f <- exp.times(t_end, f_rate)$t
  v_rate <- -log(0.5)/halflife
  sum <- 0
  for(i in 1:length(t_f)){
    ind <- which(t_c<=t_f[i])
    if(length(ind)>0){
      last <- t_c[ind[length(ind)]] 
      sum <- sum + exp(-v_rate*(t_f[i]-last))
    }
  }
  return(-log(SAR)/sum)
}


###############################
infPeriod <- 5
# Without hand washing
dc <- seq(1,30)/60 # Persistence of virus on hands (per hour)
eps <- probTrans <- meanP <- rep(0,length(dc))
for(i in 1:length(dc)){
  eps[i] <- compute.eps(t_end=24*infPeriod,c_mean=1,f_rate=11,halflife=dc[i], SAR=0.5, seed=12345)
  p <- household.fun(t_end=24*infPeriod, HW=0, f_rate=11, eps=eps[i], halflife=dc[i],seed=12345)$v
  meanP[i] <- mean(p)
}

plot(dc, eps, 
     xlab="Half life of probability \n of persistence (hours)", 
     ylab="Probability of transmission per contact",
     main="No hand washing \n SAR=0.5")

# =============================================================================== #
# how transmission prob in a 2 person household changes with hand hygiene 
# frequency as a function of duration of contamination (but where parameters are 
# constrained to give the same secondary attack rate in the absence of hand hygiene)
# =============================================================================== #
# hw <- c(1,2,3,6,12,15,20)
hw <- c(1,2,4,8,16)
dc <- seq(1,30)/60
data <- dataInf <- NULL
for(h in hw){
  eps <- probTrans <- meanP <- pInf <- rep(0,length(dc))
  for(i in 1:length(dc)){
    eps[i] <- compute.eps(t_end=24*infPeriod,c_mean=1,f_rate=11,halflife=dc[i], seed=12345)
    pInf[i] <- household.fun(t_end=24*infPeriod,HW=1,hw_mean=h, f_rate=11, eps=eps[i], halflife=dc[i],seed=12345)$p_inf
  }
  dataInf <- cbind(dataInf, pInf)
}

dataInf <- as.data.frame(dataInf)
colnames(dataInf) <- 1:ncol(dataInf)
rownames(dataInf) <- 1:nrow(dataInf)
dataInf %>% gather() %>% group_by(key) %>% 
  mutate(x=1:n()) %>%
  ggplot(aes(x=x, y=value,group=key,color=key)) + 
  geom_line() +
  xlab("Half life of probability of persistence (minutes)") + 
  ylab("Cumulative probability of infection") + 
  # scale_x_discrete(breaks=xaxis, labels=as.character(xaxis)) + 
  scale_color_discrete(name="Hand hygiene frequency \nEvery .. hours",
                       labels=hw)


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


