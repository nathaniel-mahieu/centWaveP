# centWaveP

This is a fork of the original centWave algorithm published by Tautenhan et. al. [here](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-504). The code was forked from the XCMS bioconductor package mirrored [here](https://github.com/sneumann/xcms).

Several changes have been made and this repository should be considered experimental.  Changes:

 - Separated baseline and noise estimation from the peak detection.  EICs submitted to the peak detection algorithm must include precomputed baseline and noise estimates
 - Separated ROI detection from chromatographic peak detection
 - Alternative signal to noise estimates
 - Chromatogram smoothing
 - Wavelet based valley tracking
 - Marginally improved code documentation
 - Integration bound selection is based on the supplied baseline estimates

Work is in progress. Remaining to do:

 - Modify ROI detection
 - Remove XCMS dependencies
 - Write XCMS compatibility wrapper