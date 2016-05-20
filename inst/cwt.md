


# Wavelet Scale Correspondance

```r
library(centWaveP)
```


```r
for (i in seq(1, 50, 5)) {
  wav = centWaveP:::return.wavelet('mexh') %>% { centWaveP:::scale.wavelet(i, .)$y }
  plot(wav)
  cat(sum(wav > 0), i, "\n")
  }
```

<img src="figure_cwt/unnamed-chunk-2-1.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 2 1
```

<img src="figure_cwt/unnamed-chunk-2-2.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 20 6
```

<img src="figure_cwt/unnamed-chunk-2-3.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 26 11
```

<img src="figure_cwt/unnamed-chunk-2-4.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 32 16
```

<img src="figure_cwt/unnamed-chunk-2-5.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 42 21
```

<img src="figure_cwt/unnamed-chunk-2-6.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 52 26
```

<img src="figure_cwt/unnamed-chunk-2-7.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 62 31
```

<img src="figure_cwt/unnamed-chunk-2-8.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 72 36
```

<img src="figure_cwt/unnamed-chunk-2-9.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 82 41
```

<img src="figure_cwt/unnamed-chunk-2-10.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />

```
## 92 46
```


```r
for (i in seq(1, 50, 5)) {
  wav = centWaveP:::return.wavelet('nmexh') %>% { centWaveP:::scale.wavelet(i, .)$y }
  plot(wav)
  cat(sum(wav < 0), i, "\n")
  }
```

<img src="figure_cwt/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 2 1
```

<img src="figure_cwt/unnamed-chunk-3-2.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 20 6
```

<img src="figure_cwt/unnamed-chunk-3-3.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 26 11
```

<img src="figure_cwt/unnamed-chunk-3-4.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 32 16
```

<img src="figure_cwt/unnamed-chunk-3-5.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 42 21
```

<img src="figure_cwt/unnamed-chunk-3-6.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 52 26
```

<img src="figure_cwt/unnamed-chunk-3-7.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 62 31
```

<img src="figure_cwt/unnamed-chunk-3-8.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 72 36
```

<img src="figure_cwt/unnamed-chunk-3-9.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 82 41
```

<img src="figure_cwt/unnamed-chunk-3-10.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## 92 46
```


```r
eic = { dnorm(seq(-6, 6, by =0.1)) } %>% { ./max(.) }
plot(eic)
```

<img src="figure_cwt/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```r
for (i in seq(1, 50, 5)) {
  wav = cwt(eic, i, 'mexh')
  print(plot(wav, main = i))
  cat(sum(wav > 0), i, "\n")
}
```

<img src="figure_cwt/unnamed-chunk-4-2.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 21 1
```

<img src="figure_cwt/unnamed-chunk-4-3.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 23 6
```

<img src="figure_cwt/unnamed-chunk-4-4.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 29 11
```

<img src="figure_cwt/unnamed-chunk-4-5.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 38 16
```

<img src="figure_cwt/unnamed-chunk-4-6.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 47 21
```

<img src="figure_cwt/unnamed-chunk-4-7.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 56 26
```

<img src="figure_cwt/unnamed-chunk-4-8.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 65 31
```

<img src="figure_cwt/unnamed-chunk-4-9.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 75 36
```

<img src="figure_cwt/unnamed-chunk-4-10.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 84 41
```

<img src="figure_cwt/unnamed-chunk-4-11.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## NULL
## 94 46
```


```r
eic = { dnorm(seq(-6, 6, by =0.1)) + dnorm(seq(-6, 6, by =0.1), mean = 3) } %>% { ./max(.) }
plot(eic)
```

<img src="figure_cwt/unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```r
for (i in seq(1, 50, 5)) {
  wav = cwt(eic, i, 'nmexh')
  print(plot(wav, main = i))
  cat(sum(wav < 0), i, "\n")
}
```

<img src="figure_cwt/unnamed-chunk-5-2.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 35 1
```

<img src="figure_cwt/unnamed-chunk-5-3.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 38 6
```

<img src="figure_cwt/unnamed-chunk-5-4.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 53 11
```

<img src="figure_cwt/unnamed-chunk-5-5.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 59 16
```

<img src="figure_cwt/unnamed-chunk-5-6.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 61 21
```

<img src="figure_cwt/unnamed-chunk-5-7.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 67 26
```

<img src="figure_cwt/unnamed-chunk-5-8.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 74 31
```

<img src="figure_cwt/unnamed-chunk-5-9.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 82 36
```

<img src="figure_cwt/unnamed-chunk-5-10.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 90 41
```

<img src="figure_cwt/unnamed-chunk-5-11.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

```
## NULL
## 95 46
```
