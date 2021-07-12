# ExtremeSurgeAnalysis
Statistical extreme value analysis of storm surge levels with Python.

This repository contains the code belonging to an unpublished manuscript
by Markus Reinert *et al.* (2021).  Not all the code has been uploaded
so far, but will be before the manuscript is submitted.  The reference
to the paper will be added here as soon as it is published.

Example usage is shown in the files `Time-independent_GEV_fit.py` and
`Time-dependent_GEV_fit.py`.  With the surge levels for Brest from the
GESLA-2 dataset (Woodworth et al., 2017), which need to be downloaded
separately from the website https://gesla.org/, the script
`Time-independent_GEV_fit.py` produces the following result:

![Figure of a time-independent GEV fit to extreme surge levels in Brest](results/GEV_fit_Brest.png)

Please read the header of the file `tools_GESLA.py` for more details on
how to obtain the GESLA-2 surge dataset.
