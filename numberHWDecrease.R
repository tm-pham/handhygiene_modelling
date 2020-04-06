# 
dataInf <- sar.hw(dc=5.4/60, 
                  hw=2.9/60, 
                  f_rate=10, 
                  infPeriod=0.5, 
                  c_mean=1/9, 
                  SAR.nohw=sarNoHW,
                  HW.opt=3, seed=200,it=100)
(dfReg <- dataInf$dataInf)
(num_hw_delay <- dataInf$num_hw)/12
summary(dataInf$numHW)
table(dataInf$numHW)


hc <- c(60,50,40,30,20,10,9,8,7,6,5,4,3,2,1,0.5,0.25)
hh <- c(1.15, 1.34, 1.64, 2.1, 2.97,5.1,5.54,5.9,6.48, 6.9,8.1,8.7,10.1,11.9,14,15.8,18.8)
length(hc); length(hh)
df <- data.frame(cbind(hc=hc[6:length(hh)],hh=hh[6:length(hh)]))
ggplot(data=df,aes(x=hc,y=hh,group=1)) + 
  geom_line() + 
  geom_point()

hc <- c(60,50,40,30,20,10,9,8,7,6,5,4,3,2,1)
hh_reg <- 60/c(1.15, 1.34, 1.64, 2.1, 2.97, 5.1,5.54,5.9,6.48, 6.9,8.1,8.7,10.1,11.9,14)
hh_delay <- c(0.66, 0.77, 0.93, 1.2, 1.65, 2.85,2.9,3.1,3.5,3.65,4.05,4.4,4.8,5.4,5.7)


num_delay <- rep(NA,length(hh_delay))
for(h in 1:length(hh_delay)){
  dataInf <- sar.hw(dc=10.3/60, 
                    hw=hh_delay[h]/60, 
                    f_rate=10, 
                    infPeriod=0.5, 
                    c_mean=1/hc[h], 
                    SAR.nohw=sarNoHW,
                    HW.opt=3, seed=10,it=100)
  num_delay[h] <- dataInf$num_hw/12
}



p <- plot.compare(as.data.frame(hh_reg),as.data.frame(num_delay),freq=hc, 
             col1="deepskyblue4", col2="firebrick3",
             x.title="Hand contamination rate (per hour)",
             y.title="Number of hand washes per hour \nto decrease probability of infection by 50%",
             legend.title="Hand washing",
             legend.labels=c("time-fixed","event-prompted"), point.size=2.5,points=TRUE)
p 

ggsave()