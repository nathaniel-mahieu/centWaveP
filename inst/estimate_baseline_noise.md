

# Baseline Estimations


```r
library(centWaveP)
library(ggplot2)
```

## HILIC, 60 minute gradient, 150 mm x 1 mm x 3 um, 50 uL/min

```r
eic = readRDS("eic.rds")

ggplot(eic) + 
  geom_line(aes(x = rt, y = i))
```

<img src="figure_baseline/unnamed-chunk-2-1.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```r
ggplot(subset(eic,abs(rt - 1300) < 200)) + 
  geom_line(aes(x = rt, y = i))
```

<img src="figure_baseline/unnamed-chunk-2-2.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```r
eic.n = estimateBaselineNoise(eic, peakwidth = c(15,70), minslope.peak = 2000, plot.tf = T)
```

<img src="figure_baseline/unnamed-chunk-2-3.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" /><img src="figure_baseline/unnamed-chunk-2-4.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```r
ggplot(subset(eic.n, abs(rt - 1300) < 200)) + 
  geom_line(aes(x = rt, y = i)) +
  geom_line(aes(x = rt, y =noise.local.sd), colour = "red") +
  geom_line(aes(x = rt, y =baseline), colour = "orchid") +
  geom_line(aes(x = rt, y =noise.baseline.sd), colour = "green")
```

<img src="figure_baseline/unnamed-chunk-2-5.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```r
ggplot(subset(eic.n,abs(rt - 1020) < 70)) + 
  geom_line(aes(x = rt, y = i)) +
  geom_line(aes(x = rt, y =noise.local.sd), colour = "red") +
  geom_line(aes(x = rt, y =baseline), colour = "orchid") +
  geom_line(aes(x = rt, y =noise.baseline.sd), colour = "green")
```

<img src="figure_baseline/unnamed-chunk-2-6.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

## Plateau peak

```r
eic = { dnorm(seq(-6, 6, by =0.1)) + dnorm(seq(-6, 6, by =0.1), mean = 2) } %>% { ./max(.) }
eic = data.frame( i = eic, rt = seq(eic), scan = seq(eic))

ggplot(eic) + 
  geom_line(aes(x = rt, y = i))
```

<img src="figure_baseline/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```r
eic.n = estimateBaselineNoise(eic, peakwidth = c(15, 70), minslope.peak = .005, plot.tf = T)
```

<img src="figure_baseline/unnamed-chunk-3-2.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" /><img src="figure_baseline/unnamed-chunk-3-3.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```r
ggplot(eic.n) + 
  geom_line(aes(x = rt, y = i)) +
  geom_line(aes(x = rt, y =noise.local.sd), colour = "red") +
  geom_line(aes(x = rt, y =baseline), colour = "orchid") +
  geom_line(aes(x = rt, y =noise.baseline.sd), colour = "green")
```

<img src="figure_baseline/unnamed-chunk-3-4.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```r
eic.n = estimateBaselineNoise(eic, peakwidth = c(15, 20), minslope.peak = .005, plot.tf = T)
```

<img src="figure_baseline/unnamed-chunk-3-5.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" /><img src="figure_baseline/unnamed-chunk-3-6.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```r
ggplot(eic.n) + 
  geom_line(aes(x = rt, y = i)) +
  geom_line(aes(x = rt, y =noise.local.sd), colour = "red") +
  geom_line(aes(x = rt, y =baseline), colour = "orchid") +
  geom_line(aes(x = rt, y =noise.baseline.sd), colour = "green") + main("Inappropriate max peakwidth setting.")
```

```
## Error in eval(expr, envir, enclos): could not find function "main"
```
