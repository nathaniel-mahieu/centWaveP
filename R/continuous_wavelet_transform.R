#' Apply the continuous wavelet transform to a signal
#' 
#' \code{cwt} returns a matrix of wavelet coefficients for the specified scales.
#' 
#' This function is a modification of the original cwt code provided by the centWave algorithm [1], [2] and [3].
#' 
#' @param obs vector of the time dependent signal.  In our case, intensities.
#' @param scales vector of the wavelet scales which to apply.
#' @param psi the wavelet to apply. "mexh", "nmexh" or a data.frame containing the wavelet.
#' 
#' @return A matrix, rows corresponding to wavelet scales, columns to time.
#' 
#' 
#' @section Attribution:
#' This code was modified from the originally published centWave algorithm [1].  The code was orginally distributed and obtained under the GPL2 license via the xcms software package [2]. The original algorithms depend on the wavelet analysis code included in the MassSpecWavelet package [3]. All code herein was obtained under the GPL2 license and remains under the GPL 3 license or greater.
#' 
#' @seealso 
#' \link{\code{wave}} \link[dest=https://en.wikipedia.org/wiki/Continuous_wavelet_transform]{Continuous Wavelet Transform}
#' 
#' [1] Tautenhahn, R., Böttcher, C., & Neumann, S. (2008). Highly sensitive feature detection for high resolution LC/MS. BMC bioinformatics, 9(1), 1. \link[dest=http://dx.doi.org/10.1186/1471-2105-9-504]{http://dx.doi.org/10.1186/1471-2105-9-504}
#' [2] Smith, C. A., Want, E. J., O'Maille, G., Abagyan, R., & Siuzdak, G. (2006). XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification. Analytical chemistry, 78(3), 779-787. \link[dest=10.1021/ac051437y]{10.1021/ac051437y}
#' [3] Du, P., Kibbe, W. A., & Lin, S. M. (2006). Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching. Bioinformatics, 22(17), 2059-2065. \link[dest=10.1093/bioinformatics/btl355]{10.1093/bioinformatics/btl355}
#' 
#' @export
#' 


cwt = function (obs, scales, psi = 'mexh') {
  
  psi = return.wavelet(psi)
  
  coefs = matrix(ncol = length(scales), nrow = length(obs), dimnames = list(NULL, scales))
  obs = obs %>% c(rep(mean(.[1:5]), max(scales)*10), ., rep(mean(.[(length(.) - 5):length(.)]), max(scales)*10))
  
  for (scale in scales) {
    psi.s = scale.wavelet(scale, psi)
    
    f = vector(length = length(obs)) %>% { .[1:length(psi.s$y)] = psi.s$y; . }
    
    coefs[,as.character(scale)] = convolve(obs, f) %>% { c(.[(length(.) - floor(length(psi.s$y)/2) + 1):length(.)], .[1:(length(.) - floor(length(psi.s$y)/2))]) } %>%  { . * (scale)^(-0.5) } %>% { .[(max(scales)*10 + 1):(length(.) - max(scales)*10)] }
  }
  
  coefs
}

scale.wavelet = function(scale, psi) {
  psi.s = data.frame(x = 1 + floor((0:(scale * max(psi$x)))/(scale * psi$x[2]))) 
  psi.s$y = psi$y[psi.s$x] %>% { rev(. - mean(.)) }
  
  psi.s
}


return.wavelet = function(psi) {
  
  if (psi == 'mexh') {
    mex = data.frame(x = seq(-5, 5, length = 512))
    mex$y = (2/sqrt(3) * pi^(-0.25)) * (1 - mex$x^2) * exp(-mex$x^2/2)
    mex$x = mex$x - mex$x[1] 
    mex$y = mex$y/max(mex$y)
    mex
  } else if (psi == 'nmexh') {
    mex = data.frame(x = seq(-5, 5, length = 512))
    mex$y = (2/sqrt(3) * pi^(-0.25)) * (1 - mex$x^2) * exp(-mex$x^2/2)
    mex$x = mex$x - mex$x[1] 
    mex$y = mex$y/max(mex$y)
    nmex = mex; nmex$y = -nmex$y
    nmex
  } else if (is.data.frame(psi)) {
    psi
  } else {
    stop("Unrecognized Wavelet")  
    }
  
  }