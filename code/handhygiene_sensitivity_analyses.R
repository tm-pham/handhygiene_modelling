# =============================================================================== #
# Hand hygiene covid-19 model 
# =============================================================================== #
rm(list=ls())
# Set your working directory here:
setwd("/Users/tm-pham/academia/github/handhygiene_modelling/code/")
dataPath <- "/Users/tm-pham/academia/github/handhygiene_modelling/data/"
plotPath <- "/Users/tm-pham/academia/github/handhygiene_modelling/figures/"

# Packages
source("handhygiene_packages.R")
# Load functions
source("handhygiene_functions.R")
source("handhygiene_plotting_functions.R")
# =============================================================================== #
# how transmission prob in a 2 person household changes with hand hygiene 
# frequency as a function of duration of contamination (but where parameters are 
# constrained to give the same secondary attack rate in the absence of hand hygiene)
# =============================================================================== #

start_time <- Sys.time()
infPeriod <- 0.5
sarNoHW <- 0.1
hl1 <- 5.4; hl2 <- 36.1
face_touch_rate <- 10
hand_cont_rate <- 4
# Regular fixed hand washing
dataInf <- sar.hw(dc=seq(1,60)/60, 
                  hw=c(5/60,15/60, 30/60,1,2,6), 
                  f_rate=face_touch_rate, 
                  infPeriod=infPeriod, 
                  c_mean=hand_cont_rate, 
                  SAR.nohw=sarNoHW,
                  HW.opt=1, seed=100,it=100)
dfReg <- dataInf$dataInf
(num_hw_reg <- dataInf$num_hw)/(24*infPeriod)
# Event-prompted hand washing
dataInfDel <- sar.hw(dc=seq(1,60)/60, 
                     hw=c(1/60,5/60,15/60,45/60,1,2,4), 
                     f_rate=face_touch_rate, 
                     infPeriod=0.5, 
                     c_mean=hand_cont_rate, 
                     SAR.nohw=sarNoHW,
                     HW.opt=3, seed=100,it=100)
dfDelay <- dataInfDel$dataInf
(num_hw_delay <- dataInfDel$num_hw)/(24*infPeriod)
data <- list(dfReg=dfReg, dfDelay=dfDelay, dataInf=dataInf, dataInfDel=dataInfDel, 
             sarNoHW=sarNoHW, face_touch_rate=face_touch_rate, hand_cont_rate=hand_cont_rate)
save(data, file=paste0(dataPath, "sar_hw_RegDelay_SAR", sarNoHW*100,"C", 1/hand_cont_rate,"_f", face_touch_rate, "HalfDay.RData"))

end_time <- Sys.time()
end_time-start_time


load(paste0(dataPath, "sar_hw_RegDelaySAR5C4HalfDay.RData"))
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
                     name="Hand washing \nfixed time interval") + 
  scale_y_continuous(limit=c(0.,sarNoHW+0.0001),breaks=c(seq(0,sarNoHW, by = 0.02),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2)
plotReg 

pDelay <- plot.fun(dfDelay,
                   x.title = "Half-life of probability of persistence (minutes)", 
                   y.title = "",
                   legend.title = "Delay of hand \n washing after \ncontamination", 
                   legend.labels = c("1 min","5 min", "15 min", "45 min", "1 hour","2 hours","4 hours"),
                   legend.position = "bottom", 
                   line.size=1.25)

colDelay <- c("lightblue","deepskyblue","blue4","royalblue","turquoise3","lightgreen","green4")
legendDelay <- c("1 min","5 min","15 min", "45 min", "1 hour","2 hours","4 hours")
plotDelay <- pDelay + 
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  labs(title="B") +
  theme(plot.title = element_text(hjust=-0.1, size=20, face="bold")) +
  scale_color_manual(values=colDelay, 
                     labels=legendDelay,
                     name="Hand washing delay \nafter contamination", 
                     guide=guide_legend(nrow = 2, ncol=4)) +
  scale_y_continuous(limit=c(0.,sarNoHW+0.0001), breaks=c(seq(0,sarNoHW, by = 0.02),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2)
plotDelay
# round(min(dfDelay),2)

grid.arrange(plotReg,plotDelay,ncol=2)
g <- arrangeGrob(plotReg, plotDelay, ncol=2)
ggsave(file=paste0(plotPath,"sar", sarNoHW*100,"C", 1/hand_cont_rate,"_f10_hwHalflife_regDelay_halfday.pdf"),g, height=9,width=18)