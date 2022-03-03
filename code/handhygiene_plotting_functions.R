# ============================================================================ #
# Hand hygiene plot functions
# ============================================================================ #

# ============================================================================ #
# Plot function
# ============================================================================ #
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
