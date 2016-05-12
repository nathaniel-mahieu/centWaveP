wave = function(eic, peakwidth=c(20,50), valleywidth.min = 10, sensitivity = 1, smooth = T) {
  # This is the core, wavelet based peak detection code.
  # Changes from the original centWave code: separated out baseline and noise level detection. centWave actually relies heavily on these functions and those were not described in the original paper.  The input to this function now includes baseline and noise.sd at each timepoint of the EIC.  Additional differences may have been introduced.  Some bugs were fixed and clarity was increased.
  ## Other changes. Integration works for varying scan rates.  Baseline integration uses supplied baseline.  Retention time is supplied as centroid.
  
  # eic columns: i, i.sg, inroi, baseline, noise.sd
  eic = as.data.frame(eic)
  
  # Turn peakwidths into scales
  scalerange = round(peakwidth / mean(diff(eic$rt)) * 0.5)
  scales = seq(from=scalerange[1], to=scalerange[2], by=2)
  scales = round(1.5^(scalerange[1]:scalerange[2])) %>% { c(.[.<scalerange[2]], scalerange[2]) }
  min.valley.scale = round(valleywidth.min / mean(diff(eic$rt)) * 0.5)
  
  #Wavelet Analysis
  wCoefs = { if (smooth) eic$i.sg else eic$i } %>% cwt(., scales, "mexh")
  rL = wCoefs %>% MassSpecWavelet::getLocalMaximumCWT(.) %>% MassSpecWavelet::getRidge(.)
  
  vL = { if (smooth) eic$i.sg else eic$i } %>% cwt(., scales, "nmexh") %>% MassSpecWavelet::getLocalMaximumCWT(.) %>% MassSpecWavelet::getRidge(.)
  valleys = sapply(vL, length) %>% { scales[.] > min.valley.scale } %>% which(.) %>% { sapply(vL[.], '[[', 1) }
  
  #Go through ridges and choose best scale
  peaks = lapply(rL, function(x) {

    x = data.frame(ridge = x, scale = scales[seq_along(x)])
    
    # Remove poor ridges
    x = x[x$ridge %in% which(eic$inroi==1),,drop=F]
    x = { if (smooth) eic$i.sg else eic$i } %>% { x[which(.[x$ridge] - eic$baseline[x$ridge] > 1/sensitivity * eic$noise.local.sd[x$ridge]),,drop=F] }
    if (nrow(x) == 0) { return(NULL) } # Peak not in ROI
    
    # Limit the possible scales to those which don't overlap a valley
    possible.scales = matrixStats::rowAlls(outer(x$ridge - x$scale, valleys, ">") | outer(x$ridge + x$scale, valleys, "<")) %>% which
    
    # Pick the scale which captures the greatest intensity
    best.scale.col = sapply(possible.scales, function(i) {
      x[i,] %>% sum(eic$i[(.$ridge - .$scale):(.$ridge + .$scale)])
    }) %>% which.max
    
    best.scale <-  x$scale[best.scale.col]
    if (length(best.scale) < 1) {return(NULL)}
    
    midpos <- x$ridge[best.scale.col]
    lwpos <- max(min(which(eic$inroi == 1)),midpos - best.scale)
    rwpos <- min(midpos + best.scale, max(which(eic$inroi == 1)))
    
    consec.above.baseline = max(sapply(strsplit(paste(as.numeric(eic$i[lwpos:rwpos] > eic$baseline[lwpos:rwpos]), collapse=""), "0"), nchar))
    consec.above.noise = max(sapply(strsplit(paste(as.numeric(eic$i[lwpos:rwpos] - eic$baseline[lwpos:rwpos] > eic$noise.local.sd[lwpos:rwpos]), collapse=""), "0"), nchar))
    
    height.above.wavelet.baseline = { if (smooth) eic$i.sg else eic$i } %>% { mean(.[(midpos-1):(midpos+1)]) - mean(approx(c(lwpos, rwpos),c(.[lwpos],.[rwpos]), xout = (midpos-1):(midpos+1))$y) }
    sn.wavelet.baseline = height.above.wavelet.baseline/eic$noise.local.sd[midpos]

    c(
      wavelet.scale = best.scale, 
      wavelet.location = midpos, 
      wavelet.start = lwpos, 
      wavelet.end = rwpos, 
      wavelet.consec.above.baseline = consec.above.baseline, 
      wavelet.consec.above.noise.sd=consec.above.noise,
      wavelet.height.above.wavelet.baseline = height.above.wavelet.baseline,
      waVelet.sn.wavelet.baseline = sn.wavelet.baseline
    )
  }) %>% do.call(rbind, .)
  
  if (is.null(peaks) ) {return(NULL)}
  
  
  peakinfo = lapply(1:dim(peaks)[1], function(p) {
    
    lm <- xcms:::descendMin(wCoefs[,as.character(peaks[p,"wavelet.scale"])], istart= peaks[p,"wavelet.location"]) ## find minima
    gap <- all(eic$i[lm[1]:lm[2]] == 0) ## looks like we got stuck in a gap right in the middle of the peak
    if ((lm[1]==lm[2]) || gap ) { ## fall-back
      lm =  { if (smooth) eic$i.sg else eic$i } %>% xcms:::descendMinTol(., startpos=c(peaks[p,"wavelet.start"], peaks[p,"wavelet.end"]), ceiling(min(scanwidth)/2))
    }
    
    ## narrow down peak rt boundaries by skipping things below baseline
    bl = (eic$i - eic$baseline)[lm[1]:lm[2]]
    lm.l <-  xcms:::findEqualGreaterUnsorted(bl,1E-10)
    lm.r <- xcms:::findEqualGreaterUnsorted(rev(bl),1E-10)
    lm <- lm + c(lm.l - 1, - (lm.r - 1) )
    
    scans = lm[1]:lm[2]
    
    c(
      descent.rtcentroid = sum(eic$rt[scans] * eic$i[scans]) / sum(eic$i[scans]),
      descent.rtmin = eic$rt[min(scans)],
      descent.rtmax = eic$rt[max(scans)],
      descent.maxo = max(eic$i[scans]),
      descent.into = sum(diff(eic$rt[scans])*zoo::rollmean(eic$i[scans],2)),
      descent.intb = sum(diff(eic$rt[scans])*zoo::rollmean(eic$i[scans] - eic$baseline[scans],2))
    )
  }) %>% do.call(rbind, .)
  
  overlaps = { outer(peaks[,"wavelet.location"], peaks[,"wavelet.start"], ">") & outer(peaks[,"wavelet.location"], peaks[,"wavelet.end"], "<") } %>% { diag(.) = F; . } %>% matrixStats::rowAnys(.)
  
  cols.tort = c("wavelet.location", "wavelet.start", "wavelet.end")
  peaks.rt = peaks[,cols.tort,drop=F]; peaks.rt[] = eic$rt[peaks.rt]; colnames(peaks.rt) = paste(cols.tort, ".rt", sep="")
  
  as.data.frame(cbind(peakinfo, peaks, peaks.rt)[!overlaps, ,drop=F])
}
