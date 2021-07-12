"""Handling of climate index data files.

Climate timeseries in the correct format can be found on
https://psl.noaa.gov/gcos_wgsp/Timeseries/ (checked on 20 May 2020).

Download those that you are interested in and place them in the folder
specified as ‘NOAA_DIR’.  Make sure that their names appear in
‘ClimateIndex’.

Written by Markus Reinert, May 2020/2021.
"""


import os
from enum import Enum

import numpy as np


# Path to the directory with the NOAA climate indices
NOAA_DIR = "data/NOAA/"


class ClimateIndex(Enum):
    """Selector for the available climate indices.

    The name of each ClimateIndex is a common abbreviation for it
    (possibly extended by ‘_SMOOTH’); the value of each ClimateIndex is
    a short description of it.  The name is also the filename (without
    the extension ‘.txt’) under which the data can be found in the
    directory ‘NOAA_DIR’.  The file is expected to be in the PSL format,
    see ‘read_NOAA_PSL_file’ for more information.
    """
    AMO = "Atlantic Multidecadal Oscillation"
    AMO_SMOOTH = "Atlantic Multidecadal Oscillation, smoothed"
    AO = "Arctic Oscillation"
    NAO = "North Atlantic Oscillation"


def load_NOAA_data(dataset: ClimateIndex) -> dict:
    """Load the timeseries of a climate index from NOAA.

    See the documentation of the class ‘ClimateIndex’ for more
    information on the argument ‘dataset’.

    The returned object is a dictionary with the following keys:
     - name (string, abbreviation or short description of the climate index)
     - full_name (string, full name or description of the climate index)
     - year_start (int, first year of the timeseries)
     - year_end (int, last year of the timeseries)
     - t_years (array, time values in years)
     - index (array, values of the climate index)
    """
    print("Reading NOAA data for the", dataset.value)
    path = os.path.join(NOAA_DIR, dataset.name + ".txt")
    year_start, year_end, index = read_NOAA_PSL_file(path)
    print("{:9_d} records".format(np.count_nonzero(~np.isnan(index))))
    return {
        "name": dataset.name.replace("_SMOOTH", " (smoothed)"),
        "full_name": dataset.value,
        "year_start": year_start,
        "year_end": year_end,
        "t_years": np.arange(year_start, year_end + 1, 1 / 12),
        "index": index,
    }


def read_NOAA_PSL_file(filename):
    """Read data from a file in the Monthly PSL Standard Format by NOAA.

    The method returns three objects in the following order:
     1. the first year of the timeseries,
     2. the last year of the timeseris,
     3. the monthly data.
    The first two values are ints, the data is a NumPy array of floats
    with length 12 times the number of years in the timeseries,
    including the first and last year.

    Description of the file format (taken from
    https://psl.noaa.gov/gcos_wgsp/Timeseries/standard/ on 20 May 2020):

    In general, the format for PSL monthly timeseries is as follows:

    year1 yearend
    year1 valjan valfeb valmar valapr valmay valjun valjul valaug valsep valoct valnov valdec
    year2 valjan valfeb valmar valapr valmay valjun valjul valaug valsep valoct valnov valdec
    year3 valjan valfeb valmar valapr valmay valjun valjul valaug valsep valoct valnov valdec
    .....
    yearend valjan valfeb valmar valapr valmay valjun valjul valaug valsep valoct valnov valdec
    missingvalue
    infoline1
    infoline2
    .....
    (as many infolines as needed).

    Where year1 is the first year of the timeseries and yearend is the
    last year. For each year, list the year and 12 monthly values, all
    separated by spaces. Do not skip any years. If there are years with
    no data or months with no data, use a missing value, preferably one
    outside the range of data. After the last year of data, there should
    be a line with a missing value. Even if you don't have missing
    values, please enter a number (not in the range of data). -999 is
    often good to use. Any line after the missing value line is treated
    as text so anything can go here. [...]
    """
    with open(filename) as f:
        year_start, year_end = [int(year) for year in f.readline().split()]
        n_years = year_end - year_start + 1
        data = np.zeros(n_years * 12)
        for i in range(n_years):
            values = f.readline().split()
            assert values[0] == str(year_start + i),\
                "in line {} of file {}, expected year {}, not {}.".format(
                    i + 2, filename, year_start + i, values[0]
                )
            data[i*12 : (i+1)*12] = values[1:]
        # Replace missing values by NaNs
        missing_value = float(f.readline())
        data[data == missing_value] = np.nan
    return year_start, year_end, data
