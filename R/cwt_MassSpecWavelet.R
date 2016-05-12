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