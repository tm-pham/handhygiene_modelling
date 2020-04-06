# =============================================================================== #
# Hand hygiene covid-19 model 
# =============================================================================== #
# Packages
library(ggplot2)
library(tidyr)
library(dplyr)

# ========================================================= #
# Generate times according to exponential distribution
# Exponentially distributed time-to-event until time t_end
exp.times <- function(t_end, rate, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
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

household.fun <- function(t_end=24*1, 
                          CT = 2, c_mean=1, c_sd=1, halflife=1, 
                          f_rate=10, eps=0.01,
                          HW = 1, hw_mean=10/60, t_delay=5/60, e_HW=1,
                          seed=12345){
  if(!is.null(seed)) set.seed(seed)
  

  
  #============================# 
  # Hand contamination times
  if(CT==1){ # Fixed times (constant)
    t_c<- seq(0, t_end, by = c_mean)
  }else if(CT==2){ # Exponential distrbution
    t_c <- exp.times(t_end,1/c_mean,seed)$t
    while(length(t_c)==0){
      seed <- seed+10
      t_c <- exp.times(t_end,1/c_mean,seed)$t
    }
  }else if(CT==3){ # Log-normal (heavy-tail)
    t_c <- lognormal.times(t_end,c_mean,c_sd)$t
  }#===========================# 
  
  #============================# 
  # Hand washing times/events
  t_HW <-NULL
  if(HW==1){# Hand washing at fixed times
    t_HW <- seq(0, t_end, by = hw_mean)
  }else if(HW==2){# Hand washing at random times (Poisson)
    t_HW <- exp.times(t_end, 1/hw_mean,seed+1)$t
  }else if(HW==3){# Event-prompted hand washing
    # t_HW <- t_c + t_delay
    # t_HW <- t_HW[which(t_HW<=t_end)]
    t_new <- NULL
    if(length(t_c)==1){
      t_HW <- t_c + t_delay
      if(t_HW>t_end) t_HW <- NULL
    }else{
      for(i in 1:(length(t_c)-1)){
        if(length(t_new)==0){
          temp <- t_c[i]+t_delay
          if(temp<=t_c[i+1]){
            t_HW <- c(t_HW, temp)
            t_new <- NULL
          }else t_new <- temp
        }else{
          if(t_new<=t_c[i+1]){
            t_HW <- c(t_HW, t_new)
            t_new <- NULL
          }
        }
      }
    }

  }#===========================# 
  
  #===========================# 
  # Face-touching events
  # Poisson distribution for events/exponential distribution for times
  t_f <- exp.times(t_end, f_rate,seed+2)$t
  #===========================# 
  
  
  # rate for exponential decay computed from halflife
  d_rate <- log(2)/halflife
  
  # Compute probability of persistence at face-touching events
  # Modeled with exponential decay
  v <- matrix(rep(0, length(t_f)*length(d_rate)),nrow=length(d_rate)) # prob. of persistence at face touching events
  # total time of exposure, i.e. cumulative time between contamination event and face-touching 
  # with no hand washing event: 
  t_exp <- last_tc <- list();
  
  for(d in 1:length(d_rate)){
    temp_exp <- temp_last_tc <- NULL
    for(i in 1:length(t_f)){
      ind <- which(t_c<=t_f[i])
      if(length(ind)>0){
        # Last contamination event/time
        last <- t_c[ind[length(ind)]] 
        if(HW==0){# No handwashing events
          v[d,i] <- exp(-d_rate[d]*(t_f[i]-last)) 
          temp_exp <- c(temp_exp, t_f[i]-last)
          temp_last_tc <- c(temp_last_tc, last)
        }else{
          # Was there a hand washing in between contamination and face-touching?
          ind_HW <- intersect(which(t_HW>=last), which(t_HW<t_f[i]))
          if(length(ind_HW)>0){# yes
            if(e_HW==1) v[d,i] <- 0 # perfect hand washing
            else{# imperfect handwashing
              tLast <- last
              tNew <- 0
              for(j in 1:length(ind_HW)){
                dt <- tNew + t_HW[ind_HW[j]]-tLast
                vNew <- (1-e_HW)*exp(-d_rate[d]*dt)
                tNew <- (log(d_rate[d])-log(vNew))/d_rate[d]
                tLast <- t_HW[ind_HW[j]]
              }
              v[d,i] <- exp(-d_rate[d]*(t_f[i]-tNew))
            }
          }else{
            v[d,i] <- exp(-d_rate[d]*(t_f[i]-last)) # no
            temp_exp <- c(temp_exp, t_f[i]-last)
            temp_last_tc <- c(temp_last_tc,last)
          }
        }
      }
    }
    t_exp[[d]] <- temp_exp
    last_tc[[d]] <- temp_last_tc
  }
  # Force of infection
  foi <- eps*v
  # Probability of infection
  p_inf <- unlist(apply(foi, 1,function(x) 1-exp(-sum(x))))
  
  if(length(d_rate)==1){
    v <- as.numeric(v)
    foi <- as.numeric(foi)
    p_inf <- as.numeric(p_inf)
    t_exp <- unlist(t_exp)
    last_tc <- unlist(last_tc)
  }
  
  return(list(p_inf=p_inf, v=v, t_c=t_c, t_f=t_f, t_HW=t_HW, t_exp=t_exp, last_tc=last_tc))
}

# Calculate epsilon (prob. of transmission per contact) in scenario with no HW
# Fixed SAR
compute.eps <- function(t_end=24,c_mean=1,f_rate=10,halflife=1/60, SAR=0.5, seed=12345){
  if(!is.null(seed)) set.seed(seed)
  d_rate <- log(2)/halflife
  t_c <- exp.times(t_end,1/c_mean, seed)$t
  while(length(t_c)==0){
    seed <- seed+10
    t_c <- exp.times(t_end,1/c_mean,seed)$t
  }
  t_f <- exp.times(t_end, f_rate, seed+2)$t
  v <- matrix(rep(0, length(t_f)*length(d_rate)),nrow=length(d_rate)) # prob. of persistence at face touching events

  for(d in 1:length(d_rate)){
    for(i in 1:length(t_f)){
      ind <- which(t_c<=t_f[i])
      if(length(ind)>0){
        last <- t_c[ind[length(ind)]] 
        v[d,i] <- exp(-d_rate[d]*(t_f[i]-last))
      }
    }
  }

  return(unlist(apply(v,1,function(x) -log(1-SAR)/sum(x))))
}



sar.hw <- function(dc=seq(1,60)/60, 
                   hw=c(5/60,15/60,0.5,1,2,6), 
                   f_rate=10, 
                   infPeriod=12, 
                   c_mean=1,CT=2, 
                   SAR.nohw=0.5,
                   HW.opt=1, 
                   seed=12345, it=100){
  dataInf <- matrix(rep(0,length(dc)*length(hw)),ncol=length(hw))
  num_hw <- rep(0, length(hw))
  numHW <- rep(0, it)
  for(h in 1:length(hw)){
    eps <- pInf <- rep(0,length(dc))
    pinf <- matrix(rep(0, it*length(dc)), nrow=it)
    i <- 1
    for(s in seed:(seed+it-1)){
      eps <- compute.eps(t_end=24*infPeriod,c_mean=c_mean,f_rate=f_rate,halflife=dc,SAR=SAR.nohw, seed=s)
      if(HW.opt%in%c(1,2)){
        sim <- household.fun(t_end=24*infPeriod,CT=CT,c_mean=c_mean,HW=HW.opt,hw_mean=hw[h],f_rate=f_rate,eps=eps,halflife=dc,seed=s)
        pinf[i,] <- sim$p_inf
        numHW[i] <- length(sim$t_HW)
      }else{
        sim <- household.fun(t_end=24*infPeriod,CT=CT,c_mean=c_mean,HW=HW.opt,t_delay=hw[h],f_rate=f_rate,eps=eps,halflife=dc,seed=s)
        pinf[i,] <- sim$p_inf
        # print(paste0("seed=",s))
        # print(paste0("pinf=",pinf[i]))
        numHW[i] <- length(sim$t_HW)
      }
      i <- i + 1
    }
    pInf <- unlist(apply(pinf, 2, mean))
    num_hw[h] <- mean(numHW) 
    dataInf[,h] <- pInf
  }
  dataInf <- as.data.frame(dataInf)
  colnames(dataInf) <- 1:ncol(dataInf)
  rownames(dataInf) <- 1:nrow(dataInf)
  return(list(dataInf=dataInf, num_hw=num_hw, numHW = numHW))
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



# Plot function
plot.fun <- function(df, SAR.nohw = NULL, x.title, y.title, legend.title, legend.labels, 
                     axis.title.x.size=25, axis.title.y.size=25, 
                     axis.text.x.size=22, axis.text.y.size=22,
                     legend.text.size=18, legend.title.size=20, legend.position = "right",
                     line.size=1.25){

  p <- df %>% gather() %>% group_by(key) %>% 
    mutate(x=1:n()) %>%
    ggplot(aes(x=x, y=value,group=key,color=key)) + 
    geom_line(size=line.size) +
    xlab(x.title) + 
    ylab(y.title) + 
    theme_bw()+
    theme(axis.title.x = element_text(size=axis.title.x.size),
          axis.title.y = element_text(size=axis.title.y.size),
          axis.text.x = element_text(size=axis.text.x.size),
          axis.text.y = element_text(size=axis.text.y.size),
          legend.text = element_text(size=legend.text.size),
          legend.title = element_text(size=legend.title.size),
          legend.position = legend.position)
  if(!is.null(SAR.nohw)){
    grob <- grobTree(textGrob("No Hand Washing", x=0.66,  y=0.97, hjust=0,
                              gp=gpar(col="black", fontsize=22, fontface="italic")))
    p <- p + geom_hline(yintercept = SAR.nohw, linetype="solid", color="black", size=line.size) + annotation_custom(grob) 
  }
  return(p)
}


plot.compare <- function(df1,df2,col1,col2,freq=seq(1,60), SAR.nohw = NULL, 
                         x.title="", y.title="Probability of infection", 
                         legend.title, legend.labels, 
                         axis.title.x.size=20, axis.title.y.size=20, 
                         axis.text.x.size=18, axis.text.y.size=18,
                         legend.text.size=16, legend.title.size=18, legend.position = "right",
                         line.size=1.25, point.size=1.25, points=FALSE){
  
  dR1 <- as.data.frame(cbind(freq=freq,df1))
  dD1 <- as.data.frame(cbind(freq=freq,df2))
  dDR1 <- merge(dR1, dD1, by="freq")
  dMelted1 <- reshape::melt(dDR1, id.var='freq')
  if(points){
    pc <- ggplot(dMelted1, aes(x=freq, y=value, col=variable)) + 
      geom_point(size=point.size) + 
      geom_line(size=line.size) +
      xlab(x.title) + 
      ylab(y.title) + 
      theme_bw() +
      theme(axis.title.x = element_text(size=axis.title.x.size),
            axis.title.y = element_text(size=axis.title.y.size),
            axis.text.x = element_text(size=axis.text.x.size),
            axis.text.y = element_text(size=axis.text.y.size),
            legend.text = element_text(size=legend.text.size),
            legend.title = element_text(size=legend.title.size),
            legend.position = legend.position) +
      scale_color_manual(values=c(col1,col2),
                         labels=legend.labels,
                         name=legend.title) 
  }else{
    pc <- ggplot(dMelted1, aes(x=freq, y=value, col=variable)) + 
      geom_line(size=line.size) +
      xlab(x.title) + 
      ylab(y.title) + 
      theme_bw() +
      theme(axis.title.x = element_text(size=axis.title.x.size),
            axis.title.y = element_text(size=axis.title.y.size),
            axis.text.x = element_text(size=axis.text.x.size),
            axis.text.y = element_text(size=axis.text.y.size),
            legend.text = element_text(size=legend.text.size),
            legend.title = element_text(size=legend.title.size),
            legend.position = legend.position) +
      scale_color_manual(values=c(col1,col2),
                         labels=legend.labels,
                         name=legend.title) 
  }

  return(pc)
}


num.delay <- function(hh_delay,hl=5.4,sar.noHW=0.5,infPeriod=0.5,f_rate=10,seed=10,it=100){
  num_delay <- rep(NA,length(hh_delay))
  for(h in 1:length(hh_delay)){
    dataInf <- sar.hw(dc=10.3/60, 
                      hw=hh_delay[h]/60, 
                      f_rate=10, 
                      infPeriod=infPeriod, 
                      c_mean=1/hc[h], 
                      SAR.nohw=sarNoHW,
                      HW.opt=3, seed=seed,it=it)
    num_delay[h] <- dataInf$num_hw/(24*infPeriod)
  }
  return(num_delay)
}
