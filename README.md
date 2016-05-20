# centWaveP
This code was modified from the originally published centWave algorithm [1].  The code was originally distributed and obtained under the GPL2 license via the xcms software package [2]. The original algorithms depend on the wavelet analysis code included in the MassSpecWavelet package [3]. All code herein was obtained under the GPL2 license and remains under the GPL 3 license or greater.

[1] Tautenhahn, R., Böttcher, C., & Neumann, S. (2008). Highly sensitive feature detection for high resolution LC/MS. BMC bioinformatics, 9(1), 1. [http://dx.doi.org/10.1186/1471-2105-9-504]
[2] Smith, C. A., Want, E. J., O'Maille, G., Abagyan, R., & Siuzdak, G. (2006). XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification. Analytical chemistry, 78(3), 779-787. [dest=10.1021/ac051437y]
[3] Du, P., Kibbe, W. A., & Lin, S. M. (2006). Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching. Bioinformatics, 22(17), 2059-2065. [dest=10.1093/bioinformatics/btl355]

## Modifications

Several changes have been made and this repository should be considered experimental.  Changes:

 - The algorithm accepts EICs rather than xcmsRaw objects for flexibility.
 - The peak detection algorithm has been separated from noise estimation code. (Noise estimation can be performed prior to wavelet analysis)
 - Wavelet based valley tracking to limit the aggregation of closely eluting peaks
 - Separated ROI detection from chromatographic peak detection
 - New peak quality measurements. (In our hands these provide more informative peak filters than the original code's signal to noise estimation.) See [/inst/comparison.md]
 - Chromatogram smoothing
 - Integration bound selection is based on the supplied baseline estimates
 - ROI detection by the cent portion of centWave allows for the continuation of ROIS across gaps.
 
 ## Goals
 
 Peak detection in LC/MS data, particularly non-ideal data such as HILIC is far from a solved problem.  The continued improvement of these algorithms is critical to the advance of the field. In addition to the improvements above, this implementation simplifies the usage of peak detection, allowing users to submit EICs rather than complex xcmsRaw objects and ROI lists to test peak detection. 
 
 ## Introductions
  - [/inst/comparison.md]
  - [/inst/cwt.md]
  - [/inst/estimate_baseline_noise.md]
  - [/inst/wave_centWave.md]