plotWavePeak = function(i, eic, peaks) {
  
  metadata = data.frame(peaks[i,,drop=F])
  
  range = c(metadata$wavelet.start.rt-50, metadata$wavelet.end.rt+50)
  ggplot(subset(data.frame(eic), rt < range[2] & rt > range[1])) + 
    geom_line(aes(x = rt, y = i)) +
    geom_line(aes(x = rt, y = baseline), colour = "orchid") +
    geom_line(aes(x = rt, y = i.sg), colour = "orange", alpha = 0.7) +
    geom_line(aes(x = rt, y = noise.local.sd), colour = "red", alpha = 0.7) +
    geom_vline(data = metadata, aes(xintercept = wavelet.start.rt), colour="blue", alpha = 0.8)+
    geom_vline(data = metadata, aes(xintercept = wavelet.end.rt), colour="blue", alpha = 0.8)+
    geom_vline(data = metadata, aes(xintercept = wavelet.location.rt), colour="blue", alpha = 0.8) + 
    geom_vline(data = metadata, aes(xintercept = descent.rtcentroid), colour="green", alpha = 0.8) + 
    geom_vline(data = metadata, aes(xintercept = descent.rtmin), colour="green", alpha = 0.8) +
    geom_vline(data = metadata, aes(xintercept = descent.rtmax), colour="green", alpha = 0.8) + 
    ggtitle(paste(metadata$sn.sg.above.min)) + 
    ylim(0, max(metadata$descent.maxo)*1.1)
}