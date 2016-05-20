#' Detect regions of interest (regions of m/z vs rt space with high densitites of peaks)
#' 
#' \code{cent.xr} A wrapper for xcmsRaw objects. See \code{\link{cent}}
#' 
#' @param ppm numeric the ppm mass error acceptable before a new ROI is initialized. This ppm is allowed on both sides of the mass.
#' @param prefilter numeric a two integer vector.  The first indicates the minimum number of peaks an ROI must contain to be retained.  The second indicates the minimum intensity those muct be.
#' @param maxskip integer The number of scans an ROI must not contain a peak before clossing the ROI.  \emph{This is an addition to the original algorithm. Useful for QE data.}
#' 
#' @return A list, one entry for each region of interest.  
#' 
#' 
#' @section Attribution:
#' This code was modified from the originally published centWave algorithm [1].  The code was orginally distributed and obtained under the GPL2 license via the xcms software package [2]. The original algorithms depend on the wavelet analysis code included in the MassSpecWavelet package [3]. All code herein was obtained under the GPL2 license and remains under the GPL 3 license or greater.
#' 
#' @seealso 
#' \link{\code{cent}}
#' 
#' [1] Tautenhahn, R., B?ttcher, C., & Neumann, S. (2008). Highly sensitive feature detection for high resolution LC/MS. BMC bioinformatics, 9(1), 1. \link[dest=http://dx.doi.org/10.1186/1471-2105-9-504]{http://dx.doi.org/10.1186/1471-2105-9-504}
#' [2] Smith, C. A., Want, E. J., O'Maille, G., Abagyan, R., & Siuzdak, G. (2006). XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification. Analytical chemistry, 78(3), 779-787. \link[dest=10.1021/ac051437y]{10.1021/ac051437y}
#' [3] Du, P., Kibbe, W. A., & Lin, S. M. (2006). Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching. Bioinformatics, 22(17), 2059-2065. \link[dest=10.1093/bioinformatics/btl355]{10.1093/bioinformatics/btl355}
#' 
#' @export
#' 
cent.xr = function(xr, ppm = 2, prefilter = c(0,0), maxskip = 0) {
  #SEXP findmzROI(SEXP mz, SEXP intensity, SEXP scanindex, SEXP mzrange, SEXP scanrange, SEXP lastscan, SEXP dev, SEXP minEntries, SEXP prefilter, SEXP noise, SEXP maxskip)
  .Call("findmzROI",
        as.double(xr@env$mz),
        as.double(xr@env$intensity),
        as.integer(xr@scanindex), 
        as.double(range(xr@env$mz)),
        as.integer(c(1,length(xr@scantime))), 
        as.integer(length(xr@scantime)),
        as.double(ppm * 1E-6), 
        as.integer(0), 
        as.integer(prefilter), 
        as.integer(0), 
        as.integer(maxskip),
        PACKAGE="centWaveP")
}

#' 
#' @export
#' 
roiToEic = function(roi, xr, padding = 50) {
  maxscan = length(xr@scantime)
  scrange = c(roi$scmin - padding, roi$scmax+padding) %>% {.[. < 1] = 1; .[. > maxscan] = maxscan; .}
  xcms::rawEIC(xr, mzrange = c(roi$mzmin, roi$mzmax), scanrange = scrange) %>% { names(.) = c("s", "i"); . } %>% { .$rt = xr@scantime[.$s]; .$inroi = .$s >= roi$scmin & .$s <= roi$scmax; . } %>% do.call(what = cbind) 
  }