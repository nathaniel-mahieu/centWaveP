#' Detect peaks in intensity vs time space by continuous wavelet transform
#' 
#' \code{wave} returns a list of all peak regions detected in an EIC as well as characteristics such as signal to noise estimations and integration regions.
#' 
#' This function is a modification of the original centWave algorithm [1]. The modifications are extensive - this code has not been tested to the extent which centWave has.
#' 
#' Modifications include:
#' \enumerate {
#' \item The algorithm accepts EICs rather than xcmsRaw objects for flexibility.
#' \item The peak detection algorithm has been separated from noise estimation code. (Noise estimation can be performed prior to analysis)
#' \item Valley tracking to limit the aggregation of closely eluting peaks
#' \item The ability to substitute smoothed EICs in relevant regions of the analysis
#' \item New peak quality measurements. (In our hands these provide more informative peak filters than the original code's signal to noise estimation.)
#' }
#' 
#' 
#' @param eic data.frame containing columns: i - the EIC intensities, i.sg - the smoothed EIC intensitites, rt - the retention time of the observation in seconds, inroi - 1 if the scan is part of the ROI being interrogated, baseline - the estimated baseline at each point, noise.sd - the estimated SD of the noise at each point.
#' @param peakwidth integer A two number vector containing the min and max expected peakwidth in seconds. \emph{Note, do not make the minimum peakwidth too small - a 20 second peakwidth can still detect 2 second peaks.  Smaller peakwidths limit the ability of the algorithm to distinguish noise from peakshape.  The only downside of larger peakwidths is discerning poorly resolved peaks.}
#' @param valleywidth.min integer The minimum width of a valley in seconds before the peak is split.  Stochastic fluctuations can generate multiple scan valleys which are challenging to discern.  It is suggested that this remain reasonably large unless there is very little noise.
#' @param sensitivity numeric Larger sensitivities will search for peaks in noisier regions as defined by noise.sd of the eic. This comes at the cost of additional processing time.
#' @param smooth boolean Substitutes \code{i} for \code{i.sg} for relevant calculations.

#' @return A matrix, rows corresponding to peaks. 
#' 
#' This matrix has many columns.  Columns prefixed with "wavelet." are analyses performed on the wavelet suggested peak bounds.  Columns prefixed with "descent." are analyses performed on the integration bounds (found by descending to the lowest point which can result in very long tails.) Peakshape comparisons like warpgroup may benefit from using the wavelet bounds as they will contain the bulk of the mass of the peak.
#' Useful columns include wavelet.fold.above.waveletbaseline and descent.fold.above.descentbaseline which have shown good stratification of peak quality and can be used to filter noise.
#' 
#' 
#' @section Attribution:
#' This code was modified from the originally published centWave algorithm [1].  The code was orginally distributed and obtained under the GPL2 license via the xcms software package [2]. The original algorithms depend on the wavelet analysis code included in the MassSpecWavelet package [3]. All code herein was obtained under the GPL2 license and remains under the GPL 3 license or greater.
#' 
#' @seealso 
#' \link{\code{wave}} \link{\code{estimateBaselineNoise}}
#' 
#' [1] Tautenhahn, R., B?ttcher, C., & Neumann, S. (2008). Highly sensitive feature detection for high resolution LC/MS. BMC bioinformatics, 9(1), 1. \link[dest=http://dx.doi.org/10.1186/1471-2105-9-504]{http://dx.doi.org/10.1186/1471-2105-9-504}
#' [2] Smith, C. A., Want, E. J., O'Maille, G., Abagyan, R., & Siuzdak, G. (2006). XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification. Analytical chemistry, 78(3), 779-787. \link[dest=10.1021/ac051437y]{10.1021/ac051437y}
#' [3] Du, P., Kibbe, W. A., & Lin, S. M. (2006). Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching. Bioinformatics, 22(17), 2059-2065. \link[dest=10.1093/bioinformatics/btl355]{10.1093/bioinformatics/btl355}
#' 
#' @export
#' 

wave = function(eic, peakwidth=c(20,50), valleywidth.min = 7, sensitivity = 1, smooth = T) {
  # This is the core, wavelet based peak detection code.
  # Changes from the original centWave code: separated out baseline and noise level detection. centWave actually relies heavily on these functions and those were not described in the original paper.  The input to this function now includes baseline and noise.sd at each timepoint of the EIC.  Additional differences may have been introduced.  Some bugs were fixed and clarity was increased.
  ## Other changes. Integration works for varying scan rates.  Baseline integration uses supplied baseline.  Retention time is supplied as centroid.
  
  # Turn peakwidths into scales
  scalerange = round(peakwidth / mean(diff(eic[,"rt"])) * 0.5)
  
  #Filter terrible EICs before slow wavelet stuff
  if (sum((eic[,"i"] - eic[,"baseline"]) > eic[,"noise.sd"], na.rm=T) < scalerange[1]) {return(NULL)}
  
  currscale = scalerange[1]
  scales = c()
  repeat {
    scales = c(scales, currscale)  
    currscale = (currscale)*1.15
    if (currscale > scalerange[2]) { scales = c(scales, scalerange[2]); break;}
    }
  scales = floor(scales)
  #scales = seq(from=scalerange[1], to=scalerange[2], by=2)
  
  min.valley.scale = round(valleywidth.min / mean(diff(eic[,"rt"])) * 0.5)
  
  #Wavelet Analysis
  wCoefs = { if (smooth) eic[,"i.sg"] else eic[,"i"] } %>% cwt(., scales, "mexh")
  rL = getLocalMaximumCWT(wCoefs) %>% getRidge(.)
  
  vL = { if (smooth) eic[,"i.sg"] else eic[,"i"] } %>% cwt(., scales, "nmexh") %>% getLocalMaximumCWT(.) %>% getRidge(.)
  if (length(vL) < 1) {
    valleys = 0
    } else {
    valleys = sapply(vL, length) %>% { scales[.] > min.valley.scale } %>% which(.) %>% { unlist(sapply(vL[.], '[[', 1)) }
    }
  
  #Go through ridges and choose best scale
  peaks = lapply(rL, function(x) {

    x = data.frame(ridge = x, scale = scales[seq_along(x)])
    
    # Remove poor ridges
    x = x[x$ridge %in% which(eic[,"inroi"]==1),,drop=F]
    x = { if (smooth) eic[,"i.sg"] else eic[,"i"] } %>% { x[which(.[x$ridge] - eic[x$ridge,"baseline"] > 1/sensitivity * eic[x$ridge,"noise.sd"]),,drop=F] }
    if (nrow(x) == 0) { return(NULL) } # Peak not in ROI
    
    # Limit the possible scales to those which don't overlap a valley
    possible.scales = matrixStats::rowAlls(outer(x$ridge - x$scale, valleys, ">") | outer(x$ridge + x$scale, valleys, "<")) %>% which
    
    # Pick the scale which captures the greatest intensity
    best.scale.col = sapply(possible.scales, function(i) {
      x[i,] %>% { c((.$ridge - .$scale), (.$ridge + .$scale)) } %>% {.[.<1] =1; .[. > length(eic[,"i"])] = length(eic[,"i"]); . }%>% sum(eic[,"i"][.])
    }) %>% which.max
    
    best.scale <-  x$scale[best.scale.col]
    if (length(best.scale) < 1) {return(NULL)}
    
    midpos <- x$ridge[best.scale.col]
    lwpos <- ceiling(max(min(which(eic[,"inroi"] == 1)),midpos - best.scale*1.2))
    rwpos <- floor(min(midpos + best.scale*1.2, max(which(eic[,"inroi"] == 1))))
    
    lwpos.ext = lwpos - 1
    rwpos.ext = rwpos + 1
    if (lwpos.ext < 1) {lwpos.ext = 1}
    if (rwpos.ext > length(eic[,"i"])) {rwpos.ext = length(eic[,"i"])}
    
    height.above.wavelet.baseline = eic[,"i"] %>% 
      { mean(.[(midpos-1):(midpos+1)], na.rm=T) - mean(approx(c(lwpos.ext, rwpos.ext),c(.[lwpos.ext],.[rwpos.ext]), xout = (midpos-1):(midpos+1))$y, na.rm=T) }
    baseline.mult = eic[,"i"] %>% 
      { mean(.[(midpos-1):(midpos+1)], na.rm=T) / mean(approx(c(lwpos.ext, rwpos.ext),c(.[lwpos.ext],.[rwpos.ext]), xout = (midpos-1):(midpos+1))$y, na.rm=T) }
    consec.above.noise = max(sapply(strsplit(paste(as.numeric(eic[,"i"][lwpos:rwpos] - mean(approx(c(lwpos.ext, rwpos.ext),c(eic[,"i"][lwpos.ext],eic[,"i"][rwpos.ext]), xout = lwpos:rwpos)$y, na.rm=T) > eic[,"noise.sd"][lwpos:rwpos]), collapse=""), "0"), nchar))
    
    c(
      wavelet.scale = best.scale, 
      wavelet.location = midpos, 
      wavelet.start = lwpos, 
      wavelet.end = rwpos, 
      
      wavelet.noise.sd = eic[,"noise.sd"][midpos],
      wavelet.height.above.waveletbaseline = height.above.wavelet.baseline,
      wavelet.fold.above.waveletbaseline = baseline.mult,
      wavelet.consec.above.noise.sd = consec.above.noise
    )
  }) %>% do.call(rbind, .)
  
  if (is.null(peaks) ) {return(NULL)}
  
  
  peakinfo = lapply(1:dim(peaks)[1], function(p) {
    
    lm <- descendMin(wCoefs[,as.character(peaks[p,"wavelet.scale"])], istart= peaks[p,"wavelet.location"]) ## find minima
    if ((abs(lm[1]-lm[2]) < scalerange[1]*0.5) || all(eic[,"i"][lm[1]:lm[2]] == 0) ) { 
      shrink = floor((peaks[p,"wavelet.end"] - peaks[p,"wavelet.start"])*.1)
      lm =  { if (smooth) eic[,"i.sg"] else eic[,"i"] } %>% descendMinTol(., startpos=c(peaks[p,"wavelet.start"]+shrink, peaks[p,"wavelet.end"]-shrink), ceiling(min(scalerange)))
    }
    
    ## narrow down peak rt boundaries by skipping things below baseline
    bl = (eic[,"i"] - eic[,"baseline"])[lm[1]:lm[2]]
    lm.l <-  findEqualGreaterUnsorted(bl,.05*max(eic[,"i"][c(peaks[p,"wavelet.start"]:peaks[p,"wavelet.end"])]))
    lm.r <- findEqualGreaterUnsorted(rev(bl),.05*max(eic[,"i"][c(peaks[p,"wavelet.start"]:peaks[p,"wavelet.end"])]))
    lm <- lm + c(lm.l - 1, - (lm.r - 1) )
    
    higher = which(eic[,"i"][lm] > eic[,"i"][c(peaks[p,"wavelet.start"], peaks[p,"wavelet.end"])])
    if (length(higher) > 0) { lm[higher] = c(peaks[p,"wavelet.start"], peaks[p,"wavelet.end"])[higher] }
    
    
    scans = lm[1]:lm[2]
    centroid.scan = sum(scans * eic[,"i"][scans]) / sum(eic[,"i"][scans])
    if (is.na(centroid.scan)) centroid.scan = mean(scans)
    centroid = sum(eic[,"rt"][scans] * eic[,"i"][scans]) / sum(eic[,"i"][scans])
    
    lwpos.ext = min(scans) - 1
    rwpos.ext = max(scans) + 1
    if (lwpos.ext < 1) {lwpos.ext = 1}
    if (rwpos.ext > length(eic[,"i"])) {rwpos.ext = length(eic[,"i"])}
    
    baseline.mult = { if (smooth) eic[,"i"] else eic[,"i"] } %>% 
      { mean(.[(centroid.scan-1):(centroid.scan+1)], na.rm=T) / mean(approx(c(lwpos.ext, rwpos.ext),c(.[lwpos.ext],.[rwpos.ext]), xout = (centroid.scan-1):(centroid.scan+1))$y, na.rm=T) }
    
    
    scans2 = if (length(scans) < 2) c(NA, scans) else scans
    
    c(
      descent.rtcentroid = centroid,
      descent.rtmin = eic[,"rt"][min(scans)],
      descent.rtmax = eic[,"rt"][max(scans)],
      descent.maxo = max(eic[,"i"][scans]),
      descent.into = sum(diff(eic[,"rt"][scans2])*zoo::rollmean(eic[,"i"][scans2],2)),
      descent.intb = sum(diff(eic[,"rt"][scans2])*zoo::rollmean(eic[,"i"][scans2] - eic[,"baseline"][scans2],2)),
      descent.fold.above.descentbaseline = baseline.mult
    )
  }) %>% do.call(rbind, .)
  
  overlaps = { outer(peaks[,"wavelet.location"], peaks[,"wavelet.start"], ">") & outer(peaks[,"wavelet.location"], peaks[,"wavelet.end"], "<") } %>% { diag(.) = F; . } %>% matrixStats::rowAnys(.)
  
  cols.tort = c("wavelet.location", "wavelet.start", "wavelet.end")
  peaks.rt = peaks[,cols.tort,drop=F]; peaks.rt[] = eic[,"rt"][peaks.rt]; colnames(peaks.rt) = paste(cols.tort, ".rt", sep="")
  
  as.data.frame(cbind(peakinfo, peaks, peaks.rt)[!overlaps, ,drop=F])
}
