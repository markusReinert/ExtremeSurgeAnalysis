# ExtremeSurgeAnalysis
Statistical extreme value analysis of storm surge levels with Python.

This repository contains the code belonging to an unpublished manuscript
by Markus Reinert *et al.* (2021).  The reference to the paper will be
added here as soon as it is published.

Note that the data used in the manuscript cannot be published here for
copyright reasons and needs to be obtained separately.  The surge
dataset of the authors can be obtained from the corresponding author
upon request.  The GESLA-2 surge dataset can be obtained from the
website https://gesla.org/, as explained in <tools_GESLA.py>.  The NAO
index data can be obtained from
https://psl.noaa.gov/gcos_wgsp/Timeseries/, as explained in
[tools_climate.py].

Given the datasets, the main results of the manuscript can be reproduced
with the scripts (Method_1_sliding_window_analysis.py) and
`Method_2_monthly_analysis.py`.  The analysis presented in the
Discussion can be reproduced with
`Outlook_winter_shift_in_other_stations.py`.  The parameter estimates of
the full time-dependent GEV model can be calculated with
`Time-dependent_GEV_fit_with_NAO.py`.

The functions and methods in this repository may also be used with other
datasets, for related analyses, or in similar studies.  In particular,
the script `advanced_GEV_analysis.py` may be useful, which implements
methods described in the book “An Introduction to Statistical Modeling
of Extreme Values” by Stuart Coles (2001).  Example usage of this
library is shown in the files `Time-independent_GEV_fit.py` and
`Time-dependent_GEV_fit.py`.  With the surge levels for Brest from the
GESLA-2 dataset (Woodworth et al., 2017), the script
`Time-independent_GEV_fit.py` produces the following result:

![Figure of a time-independent GEV fit to extreme surge levels in Brest](results/GEV_fit_Brest.png)
