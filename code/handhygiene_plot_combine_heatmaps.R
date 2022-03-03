# ============================================================================ #
# Combine different heatmaps
# ============================================================================ #
# Set your working directory here:
setwd("/Users/tm-pham/academia/github/handhygiene_modelling/code/")
dataPath <- "/Users/tm-pham/academia/github/handhygiene_modelling/data/"
plotPath <- "/Users/tm-pham/academia/github/handhygiene_modelling/figures/heatmap/"
fixedFolder <- "fixed/"
delayFolder <- "delay/"

# ---------------------------------------------------------------------------- #
# Fixed-time hand washing
# ---------------------------------------------------------------------------- #
prefix <- "heatmap_halflifeHandcont_fixed_SAR10_halfday_HW"
time <- c("5", "15", "30", "1", "2", "6")
unit <- c("min", "min", "min", "hour", "hours", "hours")
scenarios <-  paste0("HW fixed-time interval = ", time, " ", unit)
suffix <- ".RData"
files <- list.files(path=paste0(dataPath, fixedFolder))
sim <- vector(mode="list", length=length(files))
data_fixed <- NULL
skip_to_next <- FALSE
for(i in 1:length(time)){
  tryCatch({ # This is the "try part
    sim[[i]] <- mget(load(paste0(dataPath, fixedFolder, "/", prefix, time[i], unit[i], suffix)))
  }, 
  error=function(cond){
    message(paste("File", f, "could not be loaded."))
    message(cond)
    skip_to_next <<- TRUE
  },
  warning=function(cond){
    message("Warning!")
    message(cond)
    return(NULL)
  }, finally={}
  )
  if(skip_to_next){ next }
  temp <- expand.grid(X=seq(1,60), Y=seq(1,60,by=1))
  temp$Z <- c(as.matrix(sim[[i]]$data$dfReg*100, ncol=60))
  temp$scenario <- scenarios[i]
  data_fixed <- rbind(data_fixed, temp)
}

data_fixed$scenario <- factor(data_fixed$scenario, levels=scenarios)
(heatmap_fixed <- ggplot(data_fixed, aes(x=X, y=Y, fill=Z, group=scenario)) + 
                  facet_wrap(.~scenario) + 
                  scale_fill_viridis(discrete=FALSE) +
                  geom_tile() + 
                  labs(y="Hand contamination events per hour", x="Half life of probability of persistence (minutes)", 
                       fill = paste0("Probability of infection (%)\nBaseline = ", sarNoHW*100, "%")) +
                  theme_bw() + 
                  theme(strip.text = element_text(size = 18),
                        axis.title = element_text(size=28),
                        axis.text = element_text(size=20),
                        legend.title = element_text(size=18), 
                        legend.text = element_text(size=18)))
ggsave(heatmap_fixed, file=paste0(plotPath, "heatmap_fixed_combined_SAR", sim[[1]]$data$sarNoHW*100, ".pdf"), 
       width=16, height=9)

# ---------------------------------------------------------------------------- #
# Fixed-time hand washing
# ---------------------------------------------------------------------------- #
prefix <- "heatmap_halflifeHandcont_delay_SAR10_halfday_HW"
time <- c("1", "5", "15", "45", "1", "2", "4")
unit <- c("min","min", "min", "min", "hour", "hours", "hours")
scenarios <-  paste0("HW delay after cont = ", time, " ", unit)
suffix <- ".RData"
files <- list.files(path=paste0(dataPath, delayFolder))
sim <- vector(mode="list", length=length(files))
data_delay <- NULL
skip_to_next <- FALSE
for(i in 1:length(time)){
  tryCatch({ # This is the "try part
    sim[[i]] <- mget(load(paste0(dataPath, delayFolder, prefix, time[i], unit[i], suffix)))
  }, 
  error=function(cond){
    message(paste("File", f, "could not be loaded."))
    message(cond)
    skip_to_next <<- TRUE
  },
  warning=function(cond){
    message("Warning!")
    message(cond)
    return(NULL)
  }, finally={}
  )
  if(skip_to_next){ next }
  temp <- expand.grid(X=seq(1,60), Y=seq(1,60,by=1))
  temp$Z <- c(as.matrix(sim[[i]]$data$dfReg*100, ncol=60))
  temp$scenario <- scenarios[i]
  data_delay <- rbind(data_delay, temp)
}

data_delay$scenario <- factor(data_delay$scenario, levels=scenarios)
(heatmap_delay <- ggplot(data_delay, aes(x=X, y=Y, fill=Z, group=scenario)) + 
    facet_wrap(.~scenario) + 
    scale_fill_viridis(discrete=FALSE) +
    geom_tile() + 
    labs(y="Hand contamination events per hour", x="Half life of probability of persistence (minutes)", 
         fill = paste0("Probability of infection (%)\nBaseline = ", sarNoHW*100, "%")) +
    theme_bw() + 
    theme(strip.text = element_text(size = 18),
          axis.title = element_text(size=28),
          axis.text = element_text(size=20),
          legend.title = element_text(size=18), 
          legend.text = element_text(size=18)))
ggsave(heatmap_delay, file=paste0(plotPath, "heatmap_delay_combined_SAR", sim[[1]]$data$sarNoHW*100, ".pdf"), 
       width=16, height=12)

