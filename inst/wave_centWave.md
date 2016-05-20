


# Wavelet Based Peak Detection

```r
library(centWaveP)
library(ggplot2)
```
## HILIC, 60 minute gradient, 150 mm x 1 mm x 3 um, 50 uL/min

```r
eic = readRDS("eic.rds")
eic.n = estimateBaselineNoise(eic, peakwidth = c(15,70), minslope.peak = 10000, plot.tf = T)
```

<img src="figure_wave/unnamed-chunk-2-1.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" /><img src="figure_wave/unnamed-chunk-2-2.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```r
eic.n$inroi=T

ggplot(eic) + 
  geom_line(aes(x = rt, y = i))
```

<img src="figure_wave/unnamed-chunk-2-3.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```r
ggplot(subset(eic.n,abs(rt - 1300) < 200)) + 
  geom_line(aes(x = rt, y = i)) +
  geom_line(aes(x = rt, y =noise.local.sd), colour = "red") +
  geom_line(aes(x = rt, y =baseline), colour = "orchid") +
  geom_line(aes(x = rt, y =noise.baseline.sd), colour = "green")
```

<img src="figure_wave/unnamed-chunk-2-4.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```r
peaks = wave(eic.n, peakwidth = c(5,70), valleywidth.min = 10)
peaks = subset(peaks, descent.fold.above.descentbaseline > 0.95)

for (i in seq_along(peaks[,1])) {
  plotWavePeak(i, eic.n, peaks) %>% print
  #cat ("Press [enter] to continue")
  #line <- readline()
  
  #if (line == "x") {break;}
}
```

<img src="figure_wave/unnamed-chunk-2-5.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" /><img src="figure_wave/unnamed-chunk-2-6.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" /><img src="figure_wave/unnamed-chunk-2-7.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" /><img src="figure_wave/unnamed-chunk-2-8.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" /><img src="figure_wave/unnamed-chunk-2-9.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" /><img src="figure_wave/unnamed-chunk-2-10.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />


## Plateau peak

```r
eic = { dnorm(seq(-6, 6, by =0.1)) + dnorm(seq(-6, 6, by =0.1), mean = 3) } %>% { ./max(.) }
eic = data.frame( i = eic, rt = seq(eic), scan = seq(eic))
eic$inroi = T

ggplot(eic) + 
  geom_line(aes(x = rt, y = i))
```

<img src="figure_wave/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```r
eic.n = estimateBaselineNoise(eic, peakwidth = c(15, 70), minslope.peak = .005, plot.tf = T)
```

<img src="figure_wave/unnamed-chunk-3-2.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" /><img src="figure_wave/unnamed-chunk-3-3.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```r
ggplot(eic.n) + 
  geom_line(aes(x = rt, y = i)) +
  geom_line(aes(x = rt, y =noise.local.sd), colour = "red") +
  geom_line(aes(x = rt, y =baseline), colour = "orchid") +
  geom_line(aes(x = rt, y =noise.baseline.sd), colour = "green")
```

<img src="figure_wave/unnamed-chunk-3-4.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```r
peaks = wave(eic.n, peakwidth = c(5,70), valleywidth.min = 10)



for (i in seq_along(peaks[,1])) {
  plotWavePeak(i, eic.n, peaks) %>% print
}
```

<img src="figure_wave/unnamed-chunk-3-5.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" /><img src="figure_wave/unnamed-chunk-3-6.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />
