estimateBaselineNoise = function(eic, peakwidth, minslope.peak = 1000, plot.tf = F) {
  #Accepts an EIC with columns "i" and "rt" and adds columns "i.sg", "baseline.region", "baseline", "noise.sd"
  # This function wraps up baseline estimation and noise estimation.  Baseline estimation is performed on SG smoothed EIC traces.  Regions with a slope greater than minslope.peak are excluded from the baseline estimation.  The smoothed, nonpeak EIC is interpolated and taken as the baseline.  The original, nonpeak EIC is applied a rolling sd calculation and that is taken as the SD of the noise at each point.
  
  # min.peakwidth is used to set the windows for smoothing and noise SD estimation
  
  # eic is a data.frame or matrix with columns "i", "rt"
  eic = as.data.frame(eic)
  
  # Turn supplied peakwidth into scans
  scantimes = diff(eic$rt)
  if (diff(quantile(scantimes, c(0.05, 0.95)))/mean(scantimes) > 1) {
    warning("Scan rate varies by more than 1%.  Could cause errors in peak detection. Peak width used will be incorrect.")
  }
  scanwidth = round(peakwidth / mean(scantimes))
  
  # Compute a smoothed EIC
  sg.window = round(scanwidth[1]/3) %>% if (.%%2 == 0) .+1 else .
  eic$i.sg = signal::sgolayfilt(eic$i, p = 3, n = sg.window)
  
  # Find the SD of the residuals of the smoothing => Noise SD
  eic$noise.local.sd = zoo::rollapply(eic$i - eic$i.sg %>% {.[. < 0] = 0; .}, sg.window, "sd", fill=NA) %>% { .[is.na(.)] = 0; . }
  eic$noise.local.sd = zoo::rollmean(eic$noise.local.sd, round(sg.window/2), fill = NA)
  
  # Assign non-peak regions as regions with a slope lower than minslope.peak
  eic.prime = c(0,diff(eic$i.sg)/diff(eic$rt))
  
  baseline.points = zoo::rollmean(abs(eic.prime), sg.window, fill = 0) < minslope.peak
  subregions = nchar(strsplit(paste(abs(diff(baseline.points)), collapse=""), "1")[[1]])
  subregion.length.enough = subregions > scanwidth[2]*0.3
  
  eic$baseline.region = baseline.points & unlist(lapply(seq_along(subregions), function(i) { rep(x = subregion.length.enough[i], times = subregions[i]+1) }))
  
  # Interpolate baseline between peaks
  # When there is no region at the start or end use the first found region.
  intme.y = eic$i.sg; intme.y[1] = intme.y[min(which(eic$baseline.region))]; intme.y[length(intme.y)] = intme.y[max(which(eic$baseline.region))]
  intme.x = eic$baseline.region; intme.x[1] = T; intme.x[length(intme.x)] = T
  eic$baseline = approx(eic$rt[intme.x], intme.y[intme.x], xout = eic$rt)$y
  eic$baseline[eic$baseline < 0] = 0
  
  # Find SD of noise regions
  rolling.sd = zoo::rollapply(eic$i, width = sg.window, FUN = sd, fill= NA) %>% { .[is.na(.)] = 0; . }
  # When there is no region at the start or end use the first found region.
  intme.y = rolling.sd; intme.y[1] = intme.y[min(which(eic$baseline.region))]; intme.y[length(intme.y)] = intme.y[max(which(eic$baseline.region))]
  intme.x = eic$baseline.region; intme.x[1] = T; intme.x[length(intme.x)] = T
  eic$noise.baseline.sd = approx(eic$rt[intme.x], intme.y[intme.x], xout = eic$rt)$y
  eic$noise.baseline.sd[eic$noise.sd < 0] = 0
  
  # Plots of noise regions
  if (plot.tf) {
    plot(eic$rt[eic$baseline.region], eic$i[eic$baseline.region])
    plot(eic$rt[!eic$baseline.region], eic$i[!eic$baseline.region])
  }
  
  eic  # returns a data.frame with columns "i", "i.sg", "rt", "baseline.region", "baseline", "noise.baseline.sd", "noise.local.sd"
}