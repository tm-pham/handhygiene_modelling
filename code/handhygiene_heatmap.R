# =============================================================================== #
# Hand hygiene model
# For Figure S20 and Figure S21
# =============================================================================== #
rm(list=ls())
# Packages
library(ggplot2);library(tidyr);library(dplyr); library(gridExtra)
library(grid); library(hrbrthemes); library(viridis)

# Set your working directory here:
setwd("/Users/tm-pham/surfdrive/PHD/Oxford/covid-19/handhygiene")
dataPath <- "/Users/tm-pham/surfdrive/PHD/Oxford/covid-19/handhygiene/data/"
plotPath <- "/Users/tm-pham/surfdrive/PHD/Oxford/covid-19/handhygiene/figures/"

# Load functions
source("covid-19_handhygiene_functions.R")
# =============================================================================== #
# how transmission prob in a 2 person household changes with hand hygiene 
# frequency as a function of duration of contamination (but where parameters are 
# constrained to give the same secondary attack rate in the absence of hand hygiene)
# =============================================================================== #

start_time <- Sys.time()
infPeriod <- 0.5
sarNoHW <- 0.1
hl1 <- 5.4; hl2 <- 36.1
halflife_vec <- seq(1,60)/60
hand_cont_rate <- 1/seq(0,60,by=1)
hand_cont_rate <- hand_cont_rate[!is.infinite(hand_cont_rate)]


for(i in c(1,5,15,30)){
  handwash_fixed <- i/60
if(handwash_fixed*60<1){
  handwash_text <- handwash_fixed*60
  suffix <- "min"
}else{
  handwash_text <- handwash_fixed
  if(handwash_text<2){
    suffix <- "hour"
  }else suffix <- "hours"
}
# Regular fixed hand washing
dataInf <- sar.hw(dc=halflife_vec, 
                  hw=handwash_fixed, 
                  f_rate=10, 
                  infPeriod=infPeriod, 
                  c_mean=hand_cont_rate, 
                  SAR.nohw=sarNoHW,
                  HW.opt=1, seed=100,it=100)
dfReg <- dataInf$dataInf
(num_hw_reg <- dataInf$num_hw)/(24*infPeriod)

data <- list(dfReg=dfReg, dataInf=dataInf, 
             infPeriod=infPeriod, sarNoHW = sarNoHW, 
             halflife_vec=halflife_vec, 
             hand_cont_rate = hand_cont_rate, handwash=handwash_fixed)
save(data, file=paste0(dataPath, "heatmap_halflifeHandcont_fixed_SAR", sarNoHW*100, "_halfday_HW", handwash_text,suffix, ".RData"))

data_fixed <- expand.grid(X=seq(1,60), Y=seq(1,60,by=1))
data_fixed$Z <- c(as.matrix(dfReg*100, ncol=60))


(heatmap_fixed <- ggplot(data_fixed, aes(X, Y, fill= Z)) + 
                  geom_tile() + 
                  labs(y="Hand contamination events per hour", x="Half life of probability of persistence (minutes)", 
                  fill = paste0("Probability of infection (%)\nBaseline = ", sarNoHW*100, "%"),
                  title= paste0("Hand washing time interval (fixed) = ", handwash_text, " ",suffix)) +
                  scale_fill_viridis(discrete=FALSE) +
                  scale_x_continuous(breaks=seq(0,60, by=10)) + 
                  scale_y_continuous(breaks=seq(0,60, by=10)) +
                  theme_bw() + 
                  theme(title = element_text(size=18),
                        axis.title = element_text(size=25),
                        axis.text = element_text(size=18),
                        legend.title = element_text(size=18), 
                        legend.text = element_text(size=16)))
ggsave(heatmap_fixed, 
       file=paste0(plotPath, "heatmap/heatmap_fixed_SAR", sarNoHW*100, "_halfday_HW", handwash_text, suffix, ".pdf"), 
       width=16, height=9)
}
# Event-prompted hand washing
handwash_delay <- 4
if(handwash_delay<1){
  handwash_text <- handwash_delay*60
  suffix <- "min"
}else{
  handwash_text <- handwash_delay
  if(handwash_text<2){
    suffix <- "hour"
  }else suffix <- "hours"
}
dataInfDel <- sar.hw(dc=halflife_vec, 
                     hw=handwash_delay, 
                     f_rate=10, 
                     infPeriod=0.5, 
                     c_mean=hand_cont_rate, 
                     SAR.nohw=sarNoHW,
                     HW.opt=3, seed=100,it=100)
dfDelay <- dataInfDel$dataInf
(num_hw_delay <- dataInfDel$num_hw)/(24*infPeriod)

data_delay <- expand.grid(X=seq(1,60), Y=seq(1,60,by=1))
data_delay$Z <- c(as.matrix(dfDelay*100, ncol=60))
(heatmap_delay <- ggplot(data_delay, aes(X, Y, fill= Z)) + 
                  geom_tile() + 
                  labs(y="Hand contamination events per hour", x="Half life of probability of persistence (minutes)", 
                       fill = paste0("Probability of infection (%)\nBaseline = ", sarNoHW*100, "%"), 
                       title= paste0("Hand washing delay after contamination = ", handwash_text, " ",suffix)) +
                  scale_fill_viridis(discrete=FALSE) +
                  scale_x_continuous(breaks=seq(0,60, by=10)) + 
                  scale_y_continuous(breaks=seq(0,60, by=10)) +
                  theme_bw() + 
                  theme(title = element_text(size=18),
                        axis.title = element_text(size=25),
                        axis.text = element_text(size=18),
                        legend.title = element_text(size=18), 
                        legend.text = element_text(size=16)))
ggsave(heatmap_delay, 
       file=paste0(plotPath, "heatmap/heatmap_delay_SAR", sarNoHW*100, "_halfday_HW", handwash_text, suffix, ".pdf"), 
       width=16, height=9)

data <- list(dfReg=dfDelay, dataInf=dataInf, 
             infPeriod=infPeriod, sarNoHW = sarNoHW, 
             halflife_vec=halflife_vec, 
             hand_cont_rate = hand_cont_rate, handwash=handwash_delay)
save(data, file=paste0(dataPath, "heatmap_halflifeHandcont_delay_SAR", sarNoHW*100, "_halfday_HW", handwash_text,suffix, ".RData"))

end_time <- Sys.time()
end_time-start_time


load(paste0(dataPath, "sar_hw_RegDelaySAR50C60HalfDay.RData"))
dfReg <- data$dfReg
dfDelay <- data$dfDelay
num_hw_reg <- (data$dataInf)$num_hw
num_hw_delay <- (data$dataInfDel)$num_hw


# PLOTS
pReg <- plot.fun(dfReg,
                 x.title = "Half-life of probability of persistence (minutes)", 
                 y.title = "Probability of infection",
                 legend.title = "Hand washing \ntime interval", 
                 legend.labels = c("5 min", "15 min", "30 min", "1 hour", "2 hours", "6 hours"),
                 legend.position = "bottom",
                 line.size=1.25)
colReg <- c("darkgoldenrod1","darkorange1","coral3","darkred","lightpink","orchid","purple3","midnightblue")
legendReg <- c("5 min", "15 min", "30 min", "1 hour", "2 hours", "6 hours")
y.pos <- 0.08
text.col <- "gray21"
verticalText1 <- grobTree(textGrob(bquote("H3N2 (2"~mu*"L)"), x=0.14,  y=y.pos, hjust=0,
                                   gp=gpar(col=text.col, fontsize=18, fontface="italic")))
verticalText2 <- grobTree(textGrob(bquote("H3N2 (30"~mu*"L)"), x=0.605,  y=y.pos, hjust=0,
                                   gp=gpar(col=text.col, fontsize=18, fontface="italic")))
plotReg <- pReg + 
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  labs(title="A") + 
  theme(plot.title = element_text(hjust=-0.1, size=20, face="bold")) +
  scale_color_manual(values=colReg, labels=legendReg,
                     name="Hand washing \ntime interval") + 
  scale_y_continuous(limit=c(0.,sarNoHW+0.002),breaks=c(seq(0,sarNoHW, by = 0.05),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2)
plotReg 


pDelay <- plot.fun(dfDelay,
                   x.title = "Half-life of probability of persistence (minutes)", 
                   y.title = "",
                   legend.title = "Delay of hand \n washing after \ncontamination", 
                   legend.labels = c("1 min","5 min", "15 min", "30 min", "1 hour","2 hours","4 hours"),
                   legend.position = "bottom", 
                   line.size=1.25)

colDelay <- c("lightblue","deepskyblue","blue4","royalblue","turquoise3","lightgreen","green4")
legendDelay <- c("1 min","5 min","15 min", "30 min", "1 hour","2 hours","4 hours")
plotDelay <- pDelay + 
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  labs(title="B") +
  theme(plot.title = element_text(hjust=-0.1, size=20, face="bold")) +
  scale_color_manual(values=colDelay, 
                     labels=legendDelay,
                     name="Hand washing delay \nafter contamination", 
                     guide=guide_legend(nrow = 2, ncol=4)) +
  scale_y_continuous(limit=c(0.,sarNoHW+0.002), breaks=c(seq(0,sarNoHW, by = 0.05),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2)
plotDelay
# round(min(dfDelay),2)

grid.arrange(plotReg,plotDelay,ncol=2)
g <- arrangeGrob(plotReg, plotDelay, ncol=2)
ggsave(file=paste0(plotPath,"sar5C60_hwHalflife_regDelay_halfday.pdf"),g, height=9,width=18)


# Plot number of hand washing events per hour
freq <- legendReg
numReg <- as.data.frame(cbind(freq,count=num_hw_reg))
numReg$count <- as.numeric(as.character(num_hw_reg))/(24*infPeriod)

# numExp <- as.data.frame(cbind(freq,count=num_hw_exp))
# numExp$count <- as.numeric(as.character(num_hw_exp))/(24*infPeriod)

delay <- legendDelay
numDel <- as.data.frame(cbind(freq=delay,count=num_hw_delay))
numDel$count <- as.numeric(as.character(num_hw_delay))/(24*infPeriod)


p1 <- ggplot(data=numReg, aes(x=reorder(freq,-as.numeric(count)),y=as.numeric(count))) +
  geom_bar(stat="identity")+ 
  xlab("Hand washing frequency (fixed)") + 
  ylab("Average number of hand washes per hour") +
  geom_text(aes(label=round(count,2)), vjust=0, color="black",size=3)
# p2 <- ggplot(data=numExp, aes(x=reorder(freq,-as.numeric(count)),y=as.numeric(count))) +
#       geom_bar(stat="identity")+ 
#       xlab("Hand washing frequency (at random)") + 
#       ylab("Average number of hand washes per hour")
p3 <- ggplot(data=numDel, aes(x=reorder(freq,-as.numeric(count)),y=as.numeric(count))) +
  geom_bar(stat="identity") + 
  xlab("Delay between hand contamination and hand washing") + 
  ylab("Average number of hand washes per hour")
grid.arrange(p1,p3,ncol=2)






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


###############################
# Sensitivity analysis for different secondary attack rates
start_time <- Sys.time()
epsMat <- sens.transPerContact.SAR(dc=seq(1,30)/60,infPeriod=12,c_mean=1,f_rate=10,
                                   sar_vec=c(0.15, 0.25, 0.5, 0.75))
end_time <- Sys.time()
end_time-start_time
plot.fun(epsMat, 
         x.title = "Half-life of probability of persistence (minutes)", 
         y.title = "Probability of transmission per contact",
         legend.title = "SAR", 
         legend.labels = c("15%","25%","50%","75%"))

# Sensitivity analysis for different rates of hand contamination events
# SAR = 0.5 fixed
epsMatC <- sens.transPerContact.cont(dc=seq(1,30)/60,infPeriod=12,c_mean=1,SAR=0.5,
                                     f_rate=10, c_vec=c(1/4,1/2,1,2))
plot.fun(epsMatC, 
         x.title = "Half-life of probability of persistence (minutes)", 
         y.title = "Probability of transmission per contact",
         legend.title = "Mean time between \nhand contamination", 
         legend.labels = c("15 min","30 min","1 hour","2 hours"))

