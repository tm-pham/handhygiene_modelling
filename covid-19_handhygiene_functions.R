# =============================================================================== #
# Hand hygiene covid-19 model 
# =============================================================================== #
rm(list=ls())
# Packages
library(ggplot2)
library(tidyr)
library(dplyr)

# ========================================================= #
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
# t_end = Time period (default: 24 hours)
# c_rate = Rate at which hands of person is contaminated
# f_rate = Rate at which hands touch face (per unit of time)
# halflife = Time until probability of persistence is halved (per unit of time)
# eps = Probability of transmission per contact
# HW = Option for hand washing frequency 
# 1=fixed regular times, 2=exponential, 3=after hand contamination events with delay
# HW=1: hw_mean = regular, fixed times between hand washing
# HW=2: hw_mean = mean for exponential distribution of hand washing times
# HW=3: t_delay=delay between handwashing and contamination events
# e_HW = Efficacy of handwashing (default=1,i.e. 100% efficacious)
# CT = Option for distribution of carriage duration
# =============================================================================== #

household.fun <- function(t_end=24*12, 
                          CT = 2, c_mean=1, c_sd=1, halflife=1, 
                          f_rate=10, eps=0.01,
                          HW = 1, hw_mean=10/60, t_delay=5/60, e_HW=1,
                          seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  
  # rate for exponential decay computed from halflife
  d_rate <- log(2)/halflife
  
  #============================# 
  # Hand contamination times
  if(CT==1){ # Fixed times (constant)
    t_c<- seq(0, t_end, by = c_mean)
  }else if(CT==2){ # Exponential distrbution
    t_c <- exp.times(t_end,1/c_mean)$t
  }else if(CT==3){ # Log-normal (heavy-tail)
    t_c <- lognormal.times(t_end,c_mean,c_sd)$t
  }#===========================# 
  
  #============================# 
  # Hand washing times/events
  t_HW <-NULL
  if(HW==1){# Hand washing at fixed times
    t_HW <- seq(0, t_end, by = hw_mean)
  }else if(HW==2){# Hand washing at random times (Poisson)
    t_HW <- exp.times(t_end, 1/hw_mean)$t
  }else if(HW==3){# Event-prompted hand washing
    t_HW <- t_c + t_delay
    t_HW <- t_HW[which(t_HW<=t_end)]
  }#===========================# 
  
  #===========================# 
  # Face-touching events
  # Poisson distribution for events/exponential distribution for times
  t_f <- exp.times(t_end, f_rate)$t
  #===========================# 
  
  # Compute probability of persistence at face-touching events
  # Modeled with exponential decay
  v <- rep(0, length(t_f)) # prob. of persistence at face touching events
  # total time of exposure, i.e. cumulative time between contamination event and face-touching 
  # with no hand washing event: 
  t_exp <- last_tc <- NULL
  for(i in 1:length(t_f)){
    ind <- which(t_c<=t_f[i])
    if(length(ind)>0){
      # Last contamination event/time
      last <- t_c[ind[length(ind)]] 
      if(HW==0){# No handwashing events
        v[i] <- exp(-d_rate*(t_f[i]-last)) 
        t_exp <- c(t_exp, t_f[i]-last)
        last_tc <- c(last_tc, last)
      }else{
        # Was there a hand washing in between contamination and face-touching?
        ind_HW <- intersect(which(t_HW>=last), which(t_HW<t_f[i]))
        if(length(ind_HW)>0){# yes
          if(e_HW==1) v[i] <- 0 # perfect hand washing
          else{# imperfect handwashing
            tLast <- last
            tNew <- 0
            for(j in 1:length(ind_HW)){
              dt <- tNew + t_HW[ind_HW[j]]-tLast
              vNew <- (1-e_HW)*exp(-d_rate*dt)
              tNew <- (log(d_rate)-log(vNew))/d_rate
              tLast <- t_HW[ind_HW[j]]
            }
            v[i] <- exp(-d_rate*(t_f[i]-tNew))
          }
        }else{
          v[i] <- exp(-d_rate*(t_f[i]-last)) # no
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
compute.eps <- function(t_end,c_mean,f_rate=10,halflife, SAR=0.5, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  t_c <- exp.times(t_end,1/c_mean)$t
  t_f <- exp.times(t_end, f_rate)$t
  v <- rep(0, length(t_f)) # prob. of persistence at face touching events
  d_rate <- log(2)/halflife
  for(i in 1:length(t_f)){
    ind <- which(t_c<=t_f[i])
    if(length(ind)>0){
      last <- t_c[ind[length(ind)]] 
      v[i] <- exp(-d_rate*(t_f[i]-last))
    }
  }
  return(-log(1-SAR)/sum(v))
}



# Sensitivity analysis for different secondary attack rates
sens.transPerContact.SAR <- function(dc=seq(1,30)/60,infPeriod=12,c_mean=1,f_rate=10,
                                     sar_vec=c(0.15, 0.25, 0.5, 0.75),
                                     seed=12345){
  epsMat <- data.frame(matrix(rep(0,length(dc)*4),nrow=4,ncol=length(dc)))
  for(j in 1:length(sar_vec)){
    for(i in 1:length(dc)){
      epsMat[j,i] <- compute.eps(t_end=24*infPeriod,c_mean=1,f_rate=f_rate,halflife=dc[i], SAR=sar_vec[j], seed=seed)
    }
  }
  epsMat <- as.data.frame(t(epsMat))
  colnames(epsMat) <- 1:ncol(epsMat)
  rownames(epsMat) <- 1:nrow(epsMat)
  
  return(epsMat)
}


# Sensitivity analysis for mean times between hand contamination events
sens.transPerContact.cont <- function(dc=seq(1,30)/60,infPeriod=12,c_mean=1,SAR=0.5,
                                      f_rate=10,
                                      c_vec=c(1/4,1/2,1,2,4,8,16),
                                      seed=12345){
  epsMatC <- data.frame(matrix(rep(0,length(dc)),ncol=length(dc)))
  for(j in 1:length(c_vec)){
    for(i in 1:length(dc)){
      epsMatC[j,i] <- compute.eps(t_end=24*infPeriod,c_mean=c_vec[j],f_rate=f_rate,halflife=dc[i],SAR=SAR,seed=seed)
    }
  }
  epsMatC <- as.data.frame(t(epsMatC))
  colnames(epsMatC) <- 1:ncol(epsMatC)
  rownames(epsMatC) <- 1:nrow(epsMatC)
  
  return(epsMatC)
}


sar.hw <- function(dc=seq(1,30)/60, 
                   hw=c(5/60,15/60,0.5,1,2,4,8,16), 
                   f_rate=10, 
                   infPeriod=12, 
                   c_mean=1, 
                   SAR.nohw=0.5,
                   HW.opt=1){
  dataInf <- NULL
  for(h in hw){
    eps <- probTrans <- meanP <- pInf <- rep(0,length(dc))
    for(i in 1:length(dc)){
      eps[i] <- compute.eps(t_end=24*infPeriod,c_mean=c_mean,f_rate=f_rate,halflife=dc[i],SAR=SAR.nohw, seed=12345)
      if(HW.opt%in%c(1,2)){
        pInf[i] <- household.fun(t_end=24*infPeriod,HW=HW.opt,hw_mean=h,f_rate=f_rate,eps=eps[i],halflife=dc[i],seed=12345)$p_inf
        
      }else{
        pInf[i] <- household.fun(t_end=24*infPeriod,HW=HW.opt,t_delay=h,f_rate=f_rate,eps=eps[i],halflife=dc[i],seed=12345)$p_inf
      }
    }
    dataInf <- cbind(dataInf, pInf)
  }
  return(dataInf)
}


plot.fun <- function(df, x.title, y.title, legend.title, legend.labels, 
                     axis.title.x.size=17, axis.title.y.size=17, 
                     axis.text.x.size=15, axis.text.y.size=15,
                     legend.text.size=12, legend.title.size=13){
  p <- df %>% gather() %>% group_by(key) %>% 
    mutate(x=1:n()) %>%
    ggplot(aes(x=x, y=value,group=key,color=key)) + 
    geom_line() +
    xlab(x.title) + 
    ylab(y.title) + 
    # scale_x_discrete(breaks=xaxis, labels=as.character(xaxis)) + 
    scale_color_discrete(name=legend.title,
                         labels=legend.labels) + 
    theme_bw()+
    theme(axis.title.x = element_text(size=axis.title.x.size),
          axis.title.y = element_text(size=axis.title.y.size),
          axis.text.x = element_text(size=axis.text.x.size),
          axis.text.y = element_text(size=axis.text.y.size),
          legend.text = element_text(size=legend.text.size),
          legend.title = element_text(size=legend.title.size))
  return(p)
}
