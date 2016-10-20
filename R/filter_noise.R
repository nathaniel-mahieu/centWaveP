#' Detect artifactual shoulder peaks in Orbitrap data
#' 
#' \code{findShoulderPeaks} returns a list of masses which are shoulder peaks of a supplied mass
#' 
#' @param profxr xcmsRaw object of the profile data
#' @param p named vector containng peak details such as mass, and retention time
#' @param ppm Allowance for mirroring in the mz domain
#' @param logintrat Allowance for mirroring in the intensity domain
#' @param totintrat Maximum height relative to main peak of shoulder peaks
#'
#' @return A vector of masses
#' @export
#' 
findShoulderPeaks = function(profxr, p, ppm, logintrat = 0.7, totintrat = 0.01) {
  mzrange = c(p$mz-.5, p$mz+.5)
  scans = c(which.min(abs(profxr@scantime - p$descent.rtmin)), which.min(abs(profxr@scantime - p$descent.rtmax)))
  
  spectrum = getSpecNm(profxr, mzrange = mzrange, scanrange = scans)
  
  pks <- which(diff(sign(diff(spectrum[,"intensity"], na.pad = FALSE)), na.pad = FALSE) < 0) + 2
  mzs = spectrum[pks,"mz"]
  ints = spectrum[pks,"intensity"]
  
  targetmz.i = which.min(abs(mzs-p$mz))
  targetmz = mzs[targetmz.i]
  
  rangestart = { pks[targetmz.i]-5 } %>% { if (. < 1) 1 else . }
  targeti = max(spectrum[rangestart:(min(pks[targetmz.i]+5, length(pks))),"intensity"])
  
  dmz = ppm*p$mz / 1E6
  
  wrap = mzs - targetmz
  intapp = abs(log10(outer(ints, ints, "/"))) < logintrat
  inttot = ints/targeti < totintrat
  wrapped = abs(outer(wrap, wrap, "+")) %>% { .[lower.tri(., diag=T)] =NA; . < dmz & intapp & inttot } %>% which(arr.ind = T)
  
  unique(round(spectrum[pks[wrapped],"mz"],4))
}

#' Detect lorentzian noise around tall peaks in Orbitrap data
#' \code{findLorentzianNoise} returns the indices of peaks which fall below the noise threshold
#' 
#' @param p Named vector containing peak details 
#' @param lorentzian function accepting a mass and returning the intensity of the noise at that mass
#' @param peaks The peak table to check for noise.
#'
#' @return Indices of peaks which are below the lorentzian noise threshold
#' @export
#' 
findLorentzianNoise = function(p, lorentzian, peaks) {
  which(peaks$descent.maxo < lorentzian(peaks$mz) * p$descent.maxo & peaks$descent.maxo < p$descent.maxo*0.1)
}

#' Takes a mass and resolutionand returns a lorenzian describing the noise level in that region
#' \code{getLorentzianPeakshape} returns the indices of peaks which fall below the noise threshold
#' 
#' @param mass numeric The location of the lorentzian peak
#' @param resolution Actual (mass specific) resolution at that mass
#'
#' @return A function accepting mass and returning normalized intensity
#' @export
#' 
getLorentzianPeakshape = function(mass, resolution) {
  fwhm = mass/resolution
  function(mz) dcauchy(mz, location = mass, scale = fwhm/4, log = FALSE)/dcauchy(mass, location = mass, scale = fwhm/4, log = FALSE)
}

#' Takes a peak table and filters searches for black bar artifacts and lorentzian noise
#' \code{findNoise} returns the indices of peaks which are identified as noise
#' 
#' @param peaks The peak table to filter for noise
#' @param int The intensity at which to check around a peak for noise
#' @param resolution.calc a function which accepts mass and returns the resolution
#' @param profxr an xcmsRaw object containing the profile mode data from this experiment.
#' @param mirror.ppm the ppm window in which to detect mirror artifacts
#' @param mirror.loginrat the intensity window in which to detect mirror artifacts
#' @param mirror.totinrat the maximum intensity relative to the tall peak which can be considered noise
#'
#' @return The indices in \code{peaks} which qualify as noise
#' 
#' @export
#' 

findNoise = function(peaks, int, resolution.calc, profxr, mirror.ppm=2, mirror.logintrat=0.7, mirror.totintrat=0.1) {
  
  tallpeaks = which(peaks$descent.maxo > int)
  
  foreach(i = tallpeaks, j = icount(), .combine=c) %do% {
    cat("\r", round(j/length(tallpeaks)*100), "%      ")
    
    p = peaks[i,]

    lor = getLorentzianPeakshape(p$mz, resolution.calc(p$mz))
    
    noisepeaks = findLorentzianNoise(p, lor, peaks)
    
    shouldertf = rep(F, length(peaks$mz))
    if (is.prof(xr)) {
      shoulderpeaks = findShoulderPeaks(profxr, p, mirror.ppm, mirror.logintrat, mirror.totintrat)
      if (length(shoulderpeaks) > 0) {
        shouldertf = matrixStats::rowAnys(abs(outer(peaks$mz, shoulderpeaks, "-") / peaks$mz * 1E6) < mirror.ppm)
      } 
    } else { warning('File does not appear to be profile mode.  Shoulder peak search skipped.') }
    
    which(peaks$descent.rtcentroid < p$descent.rtmax & peaks$descent.rtcentroid > p$descent.rtmin & 
            (seq_along(peaks$mz) %in% noisepeaks | shouldertf))
    
  }
  
}


is.prof = function(xr, ppm = 1) {
  l = length(xr@env$mz)/2
  mzs = xr@env$mz[(l-100):(l+100)]
  
  all(quantile(diff(mzs)/mzs[-1] * 1E6,c(0.2, 0.8)) < ppm)
  }

getSpecNm = function(object, mzrange, scanrange) {
  scan.i = scanrange[1]:scanrange[2]
  
  scans <- list(scanrange[2]-scanrange[1])
  uniquemz <- numeric()
  for (i in seq(along = scan.i)) {
    scans[[i]] <- getScan(object, scan.i[i], mzrange)
    uniquemz <- unique(c(uniquemz, scans[[i]][,"mz"]))
  }
  uniquemz <- sort(uniquemz)
  
  intmat <- matrix(nrow = length(uniquemz), ncol = length(scan.i))
  for (i in seq(along = scan.i)) {
    scan <- getScan(object, scan.i[i], mzrange)
    if (nrow(scan) < 2) {
      intmat[,i] = 0
    } else {
      intmat[,i] <- approx(scan, xout = uniquemz)$y
    }
  }
  
  points <- cbind(mz = uniquemz, intensity = rowMeans(intmat, na.rm=T))
  
  invisible(points)
}
