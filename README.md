# IAAWT
The iterated amplitude adjusted wavelet transform (IAAWT) algorithm for generating surrogate time series.

There are three Matlab codes here:
iaawt.m is the original algorithm. It is a legacy code as it used Nick Kingsbury's original toolbox (before Mathworks implemented Nick's wavelets).
If you have Nick's toolbox then it will work.

iaawt_dualtree.m is the same as the above but using the dualtree and idualtree functions in Matlab (i.e. the Mathworks' implementation of Nick's toolbox).
multiiaawt_dualtree.m is a multivariate time-series representation of the iaawt algorithm. The key difference is one has to specify (or the code does it for you) a primary time-series. The wavelet phase differences are calculated with reference to this series. I.e. the wavelet phases of this one are randomised and then the poahse differences from this series to the others are preserved.

