#==============================================#
library(reshape)

dataPath <- "/home/thi.mui.pham/covid-19/handhygiene/modelling/code/"
plotPath <- "/home/thi.mui.pham/covid-19/handhygiene/modelling/figures/"
load(paste0(dataPath, "sar_hw_RegDelaySAR10C4HalfDay.RData"))

dfReg <- data$dfReg
dfDelay <- data$dfDelay
num_hw_reg <- (data$dataInf)$num_hw
num_hw_delay <- (data$dataInfDel)$num_hw
num_hw_reg/12;num_hw_delay/12

hl1 <- 5.4; hl2 <- 36.1
colDelay <- c("lightblue","deepskyblue","blue4","royalblue","turquoise3","lightgreen","green4")
colReg <- c("darkgoldenrod1","darkorange1","coral3","darkred","lightpink","orchid")
y.pos <- 0.92
text.col <- "gray21"
verticalText1 <- grobTree(textGrob(bquote("H3N2 (2"~mu*"L)"), x=0.14,  y=y.pos, hjust=0,
                                   gp=gpar(col=text.col, fontsize=18, fontface="italic")))
verticalText2 <- grobTree(textGrob(bquote("H3N2 (30"~mu*"L)"), x=0.605,  y=y.pos, hjust=0,
                                   gp=gpar(col=text.col, fontsize=18, fontface="italic")))


ind11 <- 2; ind12 <- 1
pc1 <- plot.compare(dfReg[,ind11],dfDelay[,ind12],colReg[ind11],colDelay[ind12],
                    legend.title="Hand washing", 
                    legend.labels=c("15 min (fixed)","1 min (delay)")) +
        geom_vline(xintercept = hl1, linetype="dashed") +
        geom_vline(xintercept = hl2, linetype="dashed") +
        annotation_custom(verticalText1) + 
        annotation_custom(verticalText2) +
        # scale_y_continuous(limit=c(round(min(dfDelay),2)-0.01,sarNoHW+0.002),breaks=c(seq(round(min(dfDelay),2)-0.01,sarNoHW, by = 0.05),sarNoHW)) +
        scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) 
pc1
ind21 <- 3; ind22 <- 3
pc2 <- plot.compare(dfReg[,ind21],dfDelay[,ind22],colReg[ind21],colDelay[ind22],
                    y.title="",
                    legend.title="Hand washing", 
                    legend.labels=c("30 min (fixed)","15 min (delay)")) +
        geom_vline(xintercept = hl1, linetype="dashed") +
        geom_vline(xintercept = hl2, linetype="dashed") +
        annotation_custom(verticalText1) + 
        annotation_custom(verticalText2) +
        # scale_y_continuous(limit=c(round(min(dfDelay),2)-0.01,sarNoHW+0.002),breaks=c(seq(round(min(dfDelay),2)-0.01,sarNoHW, by = 0.05),sarNoHW)) +
        scale_x_continuous(limit=c(0.,60),breaks=seq(0,60, by = 10)) 
pc2
ind31 <- 4; ind32 <- 5
pc3 <- plot.compare(dfReg[,ind31],dfDelay[,ind32],colReg[ind31],colDelay[ind32],
                    x.title="Half-life of probability of persistence (minutes)",
                    legend.title="Hand washing", 
                    legend.labels=c("1 hour (fixed)","1 hour (delay)")) +
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
ggsave(file=paste0(plotPath,"sar50C4_compareCurves_regDelay_halfday.pdf"),gc, height=9,width=18)


# num_hw_reg <- (data$dataInf)$num_hw
# num_hw_delay <- (data$dataInfDel)$num_hw
# numReg <- as.data.frame(cbind(freq=c(5,15,30,60,240,480),num_hw_reg/(24*0.5)))
# numDelay <- as.data.frame(cbind(freq=c(1,5,15,30,60,240),num_hw_delay/(24*0.5)))
# num <- merge(numReg,numDelay,by="freq"); colnames(num) <- c("freq","Fixed","Delayed")
# numMelted <- reshape::melt(num,id.var="freq")
# pNum <- ggplot(numMelted, aes(y=value, fill=variable, x=as.factor(freq))) +
#           geom_bar(position="dodge", stat="identity") + 
#           theme_bw() + 
#           xlab("Hand washing interval/delay (minutes)") + 
#           ylab("Average number of hand washes per hour") + 
#           theme_bw() +
#           theme(axis.title.x = element_text(size=35),
#                 axis.title.y = element_text(size=35),
#                 axis.text.x = element_text(size=25),
#                 axis.text.y = element_text(size=25),
#                 legend.text = element_text(size=27),
#                 legend.title = element_text(size=30),
#                 legend.position = "right") +
#           scale_fill_discrete(labels=c("Fixed","Delayed"),
#                               name="Hand washing") 
#   
# pNum