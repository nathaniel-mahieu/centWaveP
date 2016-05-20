#
# Code obtained from [2] under the GPL2 license.
#
# [2] Smith, C. A., Want, E. J., O'Maille, G., Abagyan, R., & Siuzdak, G. (2006). XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification. Analytical chemistry, 78(3), 779-787. \link[dest=10.1021/ac051437y]{10.1021/ac051437y}
# [3] Du, P., Kibbe, W. A., & Lin, S. M. (2006). Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching. Bioinformatics, 22(17), 2059-2065. \link[dest=10.1093/bioinformatics/btl355]{10.1093/bioinformatics/btl355}
#

descendMinTol <- function(d,startpos,maxDescOutlier) {
  l <- startpos[1]; r <- startpos[2]; outl <- 0; N <- length(d)
  ## left
  while ((l > 1) && (d[l] > 0) && outl <= maxDescOutlier) {
    if (outl > 0) vpos <- opos else vpos <- l
    if (d[l-1] > d[vpos]) outl <- outl + 1 else outl <- 0
    if (outl == 1) opos <- l
    l <- l -1
  }
  if (outl > 0) l <- l + outl
  ## right
  outl <- 0;
  while ((r < N) && (d[r] > 0) && outl <= maxDescOutlier) {
    if (outl > 0) vpos <- opos else vpos <- r
    if (d[r+1] > d[vpos]) outl <- outl + 1 else outl <- 0
    if (outl == 1) opos <- r
    r <- r + 1
  }
  if (outl > 0) r <- r - outl
  c(l,r)
}


descendMin <- function(y, istart = which.max(y)) {
  
  if (!is.double(y)) y <- as.double(y)
  unlist(.C("DescendMin",
            y,
            length(y),
            as.integer(istart-1),
            ilower = integer(1),
            iupper = integer(1),
            DUP = FALSE, PACKAGE = "centWaveP")[4:5]) + 1
}

findEqualGreaterUnsorted <- function(x, value) {
  
  if (!is.double(x)) x <- as.double(x)
  .C("FindEqualGreaterUnsorted",
     x,
     length(x),
     as.double(value),
     index = integer(1),
     DUP = FALSE, PACKAGE = "centWaveP")$index + 1
}