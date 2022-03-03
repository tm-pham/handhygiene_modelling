#==============================================#
library(reshape)

dataPath <- "/Users/tm-pham/surfdrive/PHD/Oxford/covid-19/handhygiene/data/"
plotPath <- "/Users/tm-pham/surfdrive/PHD/Oxford/covid-19/handhygiene/figures/"
load(paste0(dataPath, "sar_hw_RegDelaySAR5C4HalfDay.RData"))

dfReg <- data$dfReg
dfDelay <- data$dfDelay
num_hw_reg <- (data$dataInf)$num_hw
num_hw_delay <- (data$dataInfDel)$num_hw
num_hw_reg/12;num_hw_delay/12

# ============================================================================ #
# Compare average number of hand washes per hour for hand washing schemes
# ============================================================================ #
num_hw_reg <- (data$dataInf)$num_hw
num_hw_delay <- (data$dataInfDel)$num_hw
numReg <- as.data.frame(cbind(freq=c(5,15,30,60,120,480),Fixed=num_hw_reg/(24*0.5)))
numDelay <- as.data.frame(cbind(freq=c(1,5,15,45,60,120,240),Delay=num_hw_delay/(24*0.5)))
numMelted <- as.data.frame(rbind(reshape::melt(numReg, measure.vars=c("Fixed")), 
                                 reshape::melt(numDelay, measure.vars=c("Delay"))))
numMelted$variable <- factor(numMelted$variable, levels=c("Fixed", "Delayed"))
pNum <- ggplot(numMelted, aes(y=value, fill=variable, x=as.factor(freq))) +
  geom_bar(position="dodge", stat="identity") + 
  theme_bw() + 
  xlab("Hand washing interval/delay (minutes)") + 
  ylab("Average number of hand washes per hour") + 
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=25),
        legend.position = "right") +
  scale_fill_discrete(labels=c("Fixed-time","Event-prompted"),
                      name="Hand washing") 
pNum
ggsave(pNum, file=paste0(plotPath, "compare_num_hw_C60.pdf"), 
       width=16, height=9)

# ============================================================================ #
# Compare curves for prob. of infection for two different hand washing
# schemes (comparable average number of hand washes per hour)
# ============================================================================ #
hl1 <- 5.4; hl2 <- 36.1
colDelay <- c("lightblue","deepskyblue","blue4","royalblue","turquoise3","lightgreen","green4")
colReg <- c("darkgoldenrod1","darkorange1","coral3","darkred","lightpink","orchid")
y.pos <- 0.92
text.col <- "gray21"
verticalText1 <- grobTree(textGrob(bquote("H3N2 (2"~mu*"L)"), x=0.14,  y=y.pos, hjust=0,
                                   gp=gpar(col=text.col, fontsize=18, fontface="italic")))
verticalText2 <- grobTree(textGrob(bquote("H3N2 (30"~mu*"L)"), x=0.605,  y=y.pos, hjust=0,
                                   gp=gpar(col=text.col, fontsize=18, fontface="italic")))


ind11 <- 1; ind12 <- 2
pc1 <- plot.compare(dfReg[,ind11],dfDelay[,ind12],colReg[ind11],colDelay[ind12],
                    legend.title="Hand washing", 
                    legend.labels=c("5 min (fixed)","5 min (delay)")) +
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2) +
  # scale_y_continuous(limit=c(round(min(dfDelay),2)-0.01,sarNoHW+0.002),breaks=c(seq(round(min(dfDelay),2)-0.01,sarNoHW, by = 0.05),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) 
pc1

ind21 <- 2; ind22 <- 3
pc2 <- plot.compare(dfReg[,ind21],dfDelay[,ind22],colReg[ind21],colDelay[ind22],
                    y.title="",
                    legend.title="Hand washing", 
                    legend.labels=c("15 min (fixed)","15 min (delay)")) +
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2) +
  # scale_y_continuous(limit=c(round(min(dfDelay),2)-0.01,sarNoHW+0.002),breaks=c(seq(round(min(dfDelay),2)-0.01,sarNoHW, by = 0.05),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) 
pc2

ind31 <- 4; ind32 <- 4
pc3 <- plot.compare(dfReg[,ind31],dfDelay[,ind32],colReg[ind31],colDelay[ind32],
                    x.title="Half-life of probability of persistence (minutes)",
                    legend.title="Hand washing", 
                    legend.labels=c("60 min (fixed)","45 min (delay)")) +
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2) +
  # scale_y_continuous(limit=c(round(min(dfDelay),2)-0.01,sarNoHW+0.002),breaks=c(seq(round(min(dfDelay),2)-0.01,sarNoHW, by = 0.05),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) 
pc3

ind41 <- 5; ind42 <- 6
pc4 <- plot.compare(dfReg[,ind41],dfDelay[,ind42],colReg[ind41],colDelay[ind42],
                    x.title="Half-life of probability of persistence (minutes)",
                    y.title="",
                    legend.title="Hand washing", 
                    legend.labels=c("2 hours (fixed)","2 hours (delay)")) +
  geom_vline(xintercept = hl1, linetype="dashed") +
  geom_vline(xintercept = hl2, linetype="dashed") +
  annotation_custom(verticalText1) + 
  annotation_custom(verticalText2) +
  # scale_y_continuous(limit=c(round(min(dfDelay),2)-0.01,sarNoHW+0.002),breaks=c(seq(round(min(dfDelay),2)-0.01,sarNoHW, by = 0.05),sarNoHW)) +
  scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) 
pc4

grid.arrange(pc1,pc2,pc3,pc4,ncol=2)

gc <- arrangeGrob(pc1,pc2,pc3,pc4,ncol=2)
ggsave(file=paste0(plotPath,"sar", sarNoHW*100, "C60_f10_compareCurves_regDelay_halfday.pdf"),gc, height=9,width=18)



