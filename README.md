Trial version of a "Low Memory" mass spectrometry tools using OpenMS

The main idea of this code is to explore some ideas of how some of the OpenMS
TOPP tools could work without using too much memory. Currently they all work on
the basic principle of loading all data into memory and then applying an
algorithm on it. However, not all algorithms need access to all of the data in
memory and some may only work on one spectrum / chromatogram in isolation.

This code contains several ideas on how to circumvent the problem in a
non-intrusive manner. One idea is to cache all mzML files directly onto disc
and then get read-only random-access to the spectra and chromatograms. A
different idea is to hijack the mzML Handler and handle a spectrum /
chromatogram as soon as it is read from disc.  

Both ideas are implemented here.

The following files are currently present: 

- CommonLowMemory.h contains shared functions and objects
- LowMemGaussFilter.C contains a sample low-mem "TOPP" tool

