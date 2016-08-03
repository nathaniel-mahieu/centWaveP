
#' Plot a detected peak with various additional metadata
#' 
#' \code{plotWavePeak} plots the peak, the snoothed peak, the local noise level, baseline, as well as wavelet and descent bounds
#' 
#' @param i The index of the peak to be plotted corresponding to a row in \code{peaks}
#' @param eic The eic which the peak was detected in. The input to \code{wave}
#' @param peaks The output of \code{wave}
#' 
#' @return A ggplot2 plot grob.
#' 
#' 
#' @seealso 
#' \link{\code{wave}} \link{\code{estimateBaselineNoise}}
#' 
#' @export
#' 
plotWavePeak = function(i, eic, peaks) {
  
  metadata = data.frame(peaks[i,,drop=F])
  
  eic = as.data.frame(eic)
  
  cols = c("Local SD"="red","Wavelet Bounds"="blue","Descent Bounds"="green", "Baseline" = "orchid", "Smoothed EIC" = "orange", "EIC" = "black")
  
  theme_nate <- function(base_size = 16)
  {
    theme_grey(base_size = base_size) %+replace%
      theme(
        strip.background = element_rect(colour="grey95"), 
        strip.text = element_text(base_size*0.9, colour="grey40", face="plain"),
        
        panel.background = element_blank(),
        panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_line(colour="grey99"),
        plot.title = element_text(size = rel(1.2)),
        
        axis.title.y = element_text(vjust=1, angle=90),
        axis.title = element_text(colour="grey40"),
        axis.line = element_line(colour="grey95"),
        axis.ticks = element_line(colour="grey95"),
        
        legend.position="top",
        legend.key = element_blank(),
        legend.text = element_text(colour="grey40"),
        legend.title  = element_text(colour="grey40", face="bold")
      )
  }
  
  range = c(metadata$wavelet.start.rt-25, metadata$wavelet.end.rt+25)
  ggplot(subset(eic, rt < range[2] & rt > range[1])) + 
    geom_line(aes(x = rt, y = i, col = "EIC")) +
    geom_line(aes(x = rt, y = baseline, col = "Baseline")) +
    geom_line(aes(x = rt, y = i.sg, colour = "Smoothed EIC"), alpha = 0.7) +
    geom_line(aes(x = rt, y = noise.sd, colour = "Local SD"), alpha = 0.7) +
    geom_vline(data = metadata, aes(xintercept = wavelet.start.rt+0.03, colour="Wavelet Bounds"), alpha = 0.8)+
    geom_vline(data = metadata, aes(xintercept = wavelet.end.rt+0.06, colour="Wavelet Bounds"), alpha = 0.8)+
    geom_vline(data = metadata, aes(xintercept = wavelet.location.rt+0.09, colour="Wavelet Bounds"), alpha = 0.8) +
    geom_vline(data = metadata, aes(xintercept = descent.rtcentroid, colour="Descent Bounds"), alpha = 0.8) + 
    geom_vline(data = metadata, aes(xintercept = descent.rtmin, colour="Descent Bounds"), alpha = 0.8) +
    geom_vline(data = metadata, aes(xintercept = descent.rtmax, colour="Descent Bounds"), alpha = 0.8) + 
    scale_colour_manual(name="Colours",values=cols) +
    theme_nate() + 
    coord_cartesian(ylim=c(0,max(metadata$descent.maxo)*1.1)) + ggtitle(paste(sep=" ", "df:", round(metadata$descent.fold.above.descentbaseline,1), "wf:", round(metadata$wavelet.fold.above.waveletbaseline, 1)))
}
