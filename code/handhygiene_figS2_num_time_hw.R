# =============================================================================#
# Hand hygiene model 
# =============================================================================#
# Figure S2
# ---------------------------------------------------------------------------- #
# Plot number of hand washes that are necessary to prevent 50% of transmission
# Plot time between hand contamination 

rm(list=ls())
# Set your working directory here:
setwd("/Users/tm-pham/academia/github/handhygiene_modelling/code/")
dataPath <- "/Users/tm-pham/academia/github/handhygiene_modelling/data/"
plotPath <- "/Users/tm-pham/academia/github/handhygiene_modelling/figures/"

# Load packages
source("handhygiene_packages.R")
# Load functions
source("handhygiene_functions.R")
source("handhygiene_plotting_functions.R")

infPeriod <- 12
c_mean <- 1/4
f_rate <- 10
hl1 <- 5.4/60; hl2 <- 36.1/60
sarNoHW <- 0.1
seed <- 100

hl_vec <- c(1,5.4,10,15,20,25,36.1,40,45)/60
hw_vec <- c(2,7,8,12.5,13.5,14.5,15,16,16.2)
res_vec <- list()
time_mat <- NULL
i <- 1
for(h in hl_vec){
  eps <- compute.eps(t_end=infPeriod,c_mean=c_mean,f_rate=f_rate,halflife=h,SAR=sarNoHW, seed=seed)
  res_vec[[i]] <- household.fun(t_end=infPeriod, 
                                CT=2, 
                                c_mean=c_mean, 
                                halflife=h, 
                                hw_mean=hw_vec[i]/60, 
                                eps=eps, seed=seed)
  time_mat <- rbind(time_mat, quantile(res_vec[[i]]$t_exp*60, probs=c(0.1,0.5,0.9)))
  i <- i+1
}
colnames(time_mat) <- c("quantile.0.025", "quantile.0.5","quantile.0.975")
head(time_mat)

# ---------------------------------------------------------------------------- #
# Plot of cumulative time between hand contamination and hw 
# in order to blcok 50% of transmissions
df_num_hw <- as.data.frame(cbind(halflife=hl_vec*60, num_hw=60/hw_vec))

plot_num_hw <- ggplot(df_num_hw, aes(x=halflife,y=num_hw)) + 
  geom_point(stat="identity",size=3) + 
  theme_bw() + 
  xlab("Half-life of probability of persistence (minutes)") +
  ylab("Number of hand washes \nto block 50% of transmissions") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
plot_num_hw

ggsave(file=paste0(plotPath,"sar10C4_num_hw_block50.pdf"),
       plot_num_hw, height=9,width=14)

# Plot of cumulative time between hand contamination and hw 
# in order to blcok 50% of transmissions
time_sum <- unlist(lapply(res_vec, function(x) sum(x$t_exp)))
df_time_sum <- as.data.frame(cbind(halflife=hl_vec*60, time_sum))

plot_time_sum <- ggplot(df_time_sum, aes(x=halflife,y=time_sum)) + 
  geom_point(stat="identity",size=3) + 
  theme_bw() + 
  xlab("Half-life of probability of persistence (minutes)") +
  ylab("Time between hand contamination and hand\nwashing to block 50% of transmission (hours)") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
plot_time_sum

ggsave(file=paste0(plotPath,"sar10C4_time_sum_block50.pdf"),
       plot_time_sum, height=9,width=14)


# Combined above two plots in one plot
plot_num_time <- cowplot::plot_grid(plot_num_hw, plot_time_sum)
ggsave(plot_num_time, file=paste0(plotPath, "sar10C4_num_time_block50_comb.pdf"), 
       width=18, height=8)
ggsave(plot_num_time, file=paste0(plotPath, "FigureS2.pdf"), 
       width=18, height=8)