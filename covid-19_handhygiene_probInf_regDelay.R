# ===================================================================== #
# Hand hygiene covid-19 model 
# ===================================================================== #
# Main plots and plot of average number of HW

rm(list=ls())
# Packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)

# Set your working directory here:
dataPath <- "/home/thi.mui.pham/covid-19/handhygiene/modelling/code/"
plotPath <- "/home/thi.mui.pham/covid-19/handhygiene/modelling/figures/"

# Load functions
source(paste0(dataPath, "covid-19_handhygiene_functions.R"))

# ===================================================================== #
# Probability of infection with respect to hand hygiene frequency/delay
# as a function of half-life of prob. of persistence 
# Parameters are constrained to give the same prob. of infection
# in the absence of hand hygiene
# ===================================================================== #
start_time <- Sys.time()

infPeriod <- 12             # in hours
sarNoHW <- 0.1              # baseline prob. of infection
hl1 <- 5.4; hl2 <- 36.1     # half-life values for influenza
c_mean <- 1/4               # Mean duration between hand cont. events
f_rate <- 10                # Face-touching rate
hw_reg <- c(5/60,15/60,0.5,1,2,6)           # HW frequencies (fixed-time)
hw_delay <- c(1/60,5/60,15/60,45/60,1,2,6)  # HW delays
# hw_delay=c(1/60,5/60,15/60,1,1.25,2,3.25,6)

# Fixed-time hand washing
dataInf <- sar.hw(dc=seq(1,60)/60,
                  hw=hw_reg, 
                  f_rate=f_rate, 
                  infPeriod=infPeriod, 
                  c_mean=c_mean, 
                  SAR.nohw=sarNoHW,
                  HW.opt=1, seed=100,it=100)
dfReg <- dataInf$dataInf # data frame with prob. of infection
(num_hw_reg <- dataInf$num_hw/infPeriod) # No. of HW

# Event-prompted hand washing
dataInfDel <- sar.hw(dc=seq(1,60)/60, 
                     hw=hw_delay,
                     f_rate=f_rate, 
                     infPeriod=infPeriod, 
                     c_mean=c_mean, 
                     SAR.nohw=sarNoHW,
                     HW.opt=3, seed=100,it=100)
dfDelay <- dataInfDel$dataInf
(num_hw_delay <- dataInfDel$num_hw/infPeriod)

data <- list(dfReg=dfReg, dfDelay=dfDelay, dataInf=dataInf, dataInfDel=dataInfDel)
save(data, file=paste0(plotPath, "sar_hw_RegDelaySAR10C4HalfDay.RData"))

end_time <- Sys.time()
end_time-start_time # Total run-time



# ===================================================================== #
# MAIN PLOTS
# ===================================================================== #
load(paste0(plotPath, "sar_hw_RegDelaySAR10C4HalfDay.RData"))
dfReg <- data$dfReg
dfDelay <- data$dfDelay
num_hw_reg <- (data$dataInf)$num_hw
num_hw_delay <- (data$dataInfDel)$num_hw

# Plot for fixed-time hand washing
legendReg <- c("5 min", "15 min", "30 min", "1 hour", "2 hours", "6 hours")
# legendReg <- c("1 hour","1.5 hours", "2 hours","3 hours","3.5 hours","4 hours","6 hours")
colReg <- c("darkgoldenrod1","darkorange1","coral3","darkred","lightpink","orchid")
# colReg <- c("darkgoldenrod1","darkorange1","coral3","darkred","lightpink","orchid","purple3","midnightblue")
y.pos <- 0.08; text.col <- "gray21"
verticalText1 <- grobTree(textGrob(bquote("H3N2 (2"~mu*"L)"), x=0.14,  y=y.pos, hjust=0,
                                   gp=gpar(col=text.col, fontsize=18, fontface="italic")))
verticalText2 <- grobTree(textGrob(bquote("H3N2 (30"~mu*"L)"), x=0.605,  y=y.pos, hjust=0,
                                   gp=gpar(col=text.col, fontsize=18, fontface="italic")))

pReg <- plot.fun(dfReg,
                 x.title = "Half-life of probability of persistence (minutes)", 
                 y.title = "Probability of infection",
                 legend.title = "Hand washing \ntime interval", 
                 legend.position = "bottom",
                 line.size=1.25)
plotReg <- pReg + 
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  labs(title="A") + 
  theme(plot.title = element_text(hjust=-0.1, size=20, face="bold")) +
  scale_color_manual(values=colReg, labels=legendReg,
                     name="Hand washing \ntime interval") + 
  scale_y_continuous(limit=c(0,sarNoHW),breaks=c(seq(0,sarNoHW, by = 0.02),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2)
plotReg 


legendDelay <- c("1 min","5 min","15 min", "45 min", "1 hour","2 hours","4 hours")
colDelay <- c("lightblue","deepskyblue","blue4","royalblue","turquoise3","lightgreen","green4")
# colDelay <- c("lightblue","deepskyblue","blue4","royalblue","turquoise3","lightgreen","green4","darkgreen")
pDelay <- plot.fun(dfDelay,
                   x.title = "Half-life of probability of persistence (minutes)", 
                   y.title = "",
                   legend.title = "Delay of hand \n washing after \ncontamination", 
                   legend.labels = c("1 min","5 min", "15 min", "30 min", "1 hour","2 hours","4 hours"),
                   legend.position = "bottom", 
                   line.size=1.25)
plotDelay <- pDelay + 
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  labs(title="B") +
  theme(plot.title = element_text(hjust=-0.1, size=20, face="bold")) +
  scale_color_manual(values=colDelay, 
                     labels=legendDelay,
                     name="Hand washing delay \nafter contamination", 
                     guide=guide_legend(nrow = 2, ncol=4)) +
  scale_y_continuous(limit=c(0.,sarNoHW+0.002), breaks=c(seq(0,sarNoHW, by = 0.02),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2)
plotDelay


grid.arrange(plotReg,plotDelay,ncol=2)
g <- arrangeGrob(plotReg, plotDelay, ncol=2)
ggsave(file=paste0(plotPath,"sar10C4_hwHalflife_regDelay_halfday.pdf"),g, height=9,width=18)
# ggsave(file=paste0(plotPath,"C60_Pandejpong2012_hwHalflife_regDelay_halfday.pdf"),g, height=9,width=18)
# ========END MAIN PLOTS=============================================== #

# ===================================================================== #
# Plot number of hand washing events per hour
# ===================================================================== #
numReg <- as.data.frame(cbind(freq=hw_reg,count=num_hw_reg))
numDel <- as.data.frame(cbind(freq=hw_delay,count=num_hw_delay))


p_numReg <- ggplot(data=numReg, aes(x=reorder(freq,-as.numeric(count)),y=as.numeric(count))) +
            geom_bar(stat="identity")+ 
            xlab("Hand washing frequency (fixed)") + 
            ylab("Average number of hand washes per hour") +
            geom_text(aes(label=round(count,2)), vjust=-0.5, color="black",size=7) + 
            theme_bw() + 
            theme(axis.title = element_text(size=25),
                  axis.text = element_text(size=22)) + 
            scale_x_discrete(labels=legendReg)
p_numDelay <- ggplot(data=numDel, aes(x=reorder(freq,-as.numeric(count)),y=as.numeric(count))) +
              geom_bar(stat="identity") + 
              xlab("Delay between hand contamination and hand washing") +
              geom_text(aes(label=round(count,2)), vjust=-0.5, color="black",size=7) +
              theme_bw() + 
              theme(axis.title.y = element_blank(), 
                    axis.title.x = element_text(size=25),
                    axis.text = element_text(size=22)) + 
              scale_x_discrete(labels=legendDelay)

grid.arrange(p_numReg,p_numDelay,ncol=2)
g_num <- arrangeGrob(p_numReg,p_numDelay, ncol=2)
ggsave(file=paste0(plotPath,"sar10C4_num_hw_regDelay.pdf"),
            g_num, height=9,width=18)

