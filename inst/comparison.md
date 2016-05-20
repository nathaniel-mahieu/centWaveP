

# xcms::findPeaks.centWave vs centWaveP


```r
library(centWaveP)
library(xcms)
library(gridExtra)
```

## centWaveP

```r
file = "../../../Projects/Datasets/2nm118a_nh2_negative_ecoli_1213/2NM118A_HILICAnnoWgN14_2NM111G_12.mzXML"

xr = xcmsRaw(file, profstep=0)

roi.l = cent.xr(xr, ppm = 2, prefilter = c(4,0), maxskip = 5) # ppm is added to both sides
```

```
##  % finished: 0 10 20 30 40 50 60 70 80 90 100 
##  252768 m/z ROI's.
```

```r
eic.l = lapply(sample(roi.l, 3000), roiToEic, xr, padding=100)
eic.noise.l = lapply(eic.l, estimateBaselineNoise, peakwidth = c(15, 60), minslope.peak = 20000)

peaks = lapply(seq_along(eic.noise.l), function(i) {
  wave(eic.noise.l[[i]], peakwidth = c(15, 70), valleywidth.min = 15, sensitivity = 1, smooth = T)
})

index = lapply(seq(peaks), function(i) {
  wns = peaks[[i]]$descent.fold.above.descentbaseline
  cbind(roi = rep(i, length(wns)), wsn = wns, num = seq_along(wns))
  }) %>% do.call(what = rbind)

nrow(index)
```

```
## [1] 139
```

### SN < 1

```r
ps0 = lapply(sample(which(index[,"wsn"] < 1), 16, replace = T), function(x) {
  plotWavePeak(unname(index[x, "num"]), eic.noise.l[[index[x,"roi"]]], peaks[[index[x,"roi"]]])
})
do.call(grid.arrange, ps0)  
```

<img src="figure_comparison/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />
### 1 < SN < 2

```r
ps1 = lapply(sample(index[,"wsn"] %>% { which(. > 1  & . < 2) }, 16, replace = T), function(x) {
  plotWavePeak(unname(index[x, "num"]), eic.noise.l[[index[x,"roi"]]], peaks[[index[x,"roi"]]])
})
do.call(grid.arrange, ps1)  
```

<img src="figure_comparison/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />
### 2 < SN < 3

```r
ps2 = lapply(sample(index[,"wsn"] %>% { which(. > 2  & . < 3) }, 16, replace = T), function(x) {
  plotWavePeak(unname(index[x, "num"]), eic.noise.l[[index[x,"roi"]]], peaks[[index[x,"roi"]]])
})
do.call(grid.arrange, ps2)
```

<img src="figure_comparison/unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />
### 3 < SN < 4

```r
ps3 = lapply(sample(index[,"wsn"] %>% { which(. > 3  & . < 5) }, 16, replace = T), function(x) {
  plotWavePeak(unname(index[x, "num"]), eic.noise.l[[index[x,"roi"]]], peaks[[index[x,"roi"]]])
})
do.call(grid.arrange, ps3)
```

<img src="figure_comparison/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" style="display: block; margin: auto;" />
### 5 < SN 

```r
ps4 = lapply(sample(index[,"wsn"] %>% { which(. > 5) }, 16, replace=T), function(x) {
  plotWavePeak(unname(index[x, "num"]), eic.noise.l[[index[x,"roi"]]], peaks[[index[x,"roi"]]])
})
do.call(grid.arrange, ps4)
```

<img src="figure_comparison/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" style="display: block; margin: auto;" />

## xcms::findPeaks.centWave

```r
ROI.list = xcms:::findmzROI(xr, dev = 2*1E-6, minCentroids = 4, prefilter = c(4,0), noise = 0)
```

```
##  % finished: 0 10 20 30 40 50 60 70 80 90 100 
##  319716 m/z ROI's.
```

```r
peaks.xcms = findPeaks.centWave(xr, ppm = 2, peakwidth = c(15, 60), snthresh = 0, prefilter = c(4,0), noise = 0, ROI.list = sample(ROI.list, 3000))
```

```
## 
##  Detecting chromatographic peaks ... 
##  % finished: 0 10 20 30 40 50 60 70 80 90 100 
##  530  Peaks.
```

```r
nrow(peaks.xcms)
```

```
## [1] 530
```

```r
head(peaks.xcms)
```

```
##            mz    mzmin    mzmax       rt     rtmin    rtmax       into
## [1,] 225.0079 225.0075 225.0080  115.720   92.7069  136.562 32123174.6
## [2,] 567.3541 567.3531 567.3548 1005.980  993.0220 1053.370  5647315.1
## [3,] 647.2426 647.2417 647.2433 1720.410 1710.8600 1727.740  1380589.9
## [4,] 261.0569 261.0568 261.0571  769.188  746.0230  780.455  5133386.7
## [5,] 418.1583 418.1578 418.1587 1557.680 1541.8300 1576.330  1781535.5
## [6,] 358.1093 358.1089 358.1099 1337.840 1332.7500 1345.760   205826.5
##            intb        maxo     sn
## [1,] 30641510.0 15777098.00     86
## [2,]  1415153.0   169545.69      2
## [3,]  1380574.7   192009.64 192009
## [4,]  3728750.8   257335.88      4
## [5,]  1680869.2   161874.88     14
## [6,]   204360.3    44780.32     25
```


### 1 < SN < 5

```r
 ps4 = lapply(sample(peaks.xcms[,"sn"] %>% { which(. > 1 & . < 5) }, 16, replace=T), function(i) {
   scrange = c(peaks.xcms[i,"rtmin"]-30, peaks.xcms[i,"rtmax"]+30)
   scrange[scrange < 1] = 1; scrange[scrange > max(xr@scantime)] = max(xr@scantime)
    eic = rawEIC(xr, mzrange = c(peaks.xcms[i,"mzmin"], peaks.xcms[i,"mzmax"]), rtrange = scrange)
   
    df = data.frame(eic)
    
    ggplot(df) +
      geom_line(aes(x = scan, y = intensity))
   })
do.call(grid.arrange, ps4)
```

<img src="figure_comparison/unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" style="display: block; margin: auto;" />
### 5 < SN < 10

```r
 ps4 = lapply(sample(peaks.xcms[,"sn"] %>% { which(. > 5 & . < 10) }, 16, replace=T), function(i) {
   scrange = c(peaks.xcms[i,"rtmin"]-30, peaks.xcms[i,"rtmax"]+30)
   scrange[scrange < 1] = 1; scrange[scrange > max(xr@scantime)] = max(xr@scantime)
    eic = rawEIC(xr, mzrange = c(peaks.xcms[i,"mzmin"], peaks.xcms[i,"mzmax"]), rtrange = scrange)
   
    df = data.frame(eic)
    
    ggplot(df) +
      geom_line(aes(x = scan, y = intensity))
   })
do.call(grid.arrange, ps4)
```

<img src="figure_comparison/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" style="display: block; margin: auto;" />
### 10 < SN < 100

```r
 ps4 = lapply(sample(peaks.xcms[,"sn"] %>% { which(. > 10 & . < 100) }, 16, replace=T), function(i) {
   scrange = c(peaks.xcms[i,"rtmin"]-30, peaks.xcms[i,"rtmax"]+30)
   scrange[scrange < 1] = 1; scrange[scrange > max(xr@scantime)] = max(xr@scantime)
    eic = rawEIC(xr, mzrange = c(peaks.xcms[i,"mzmin"], peaks.xcms[i,"mzmax"]), rtrange = scrange)
   
    df = data.frame(eic)
    
    ggplot(df) +
      geom_line(aes(x = scan, y = intensity))
   })
do.call(grid.arrange, ps4)
```

<img src="figure_comparison/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" style="display: block; margin: auto;" />
### SN > 100

```r
 ps4 = lapply(sample(peaks.xcms[,"sn"] %>% { which(. > 100) }, 16, replace=T), function(i) {
   scrange = c(peaks.xcms[i,"rtmin"]-30, peaks.xcms[i,"rtmax"]+30)
   scrange[scrange < 1] = 1; scrange[scrange > max(xr@scantime)] = max(xr@scantime)
    eic = rawEIC(xr, mzrange = c(peaks.xcms[i,"mzmin"], peaks.xcms[i,"mzmax"]), rtrange = scrange)
   
    df = data.frame(eic)
    
    ggplot(df) +
      geom_line(aes(x = scan, y = intensity))
   })
do.call(grid.arrange, ps4)
```

<img src="figure_comparison/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" style="display: block; margin: auto;" />
