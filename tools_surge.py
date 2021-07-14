"""Handling of skew surge data.

Two skew surge datasets are implemented: the dataset for Brest created
by Reinert et al. (2021), which can be received from the authors upon
request, and the GESLA-2 surge dataset of Marcos & Woodworth (2017).
See the file ‘tools_GESLA’ for further information on the latter.

Written by Markus Reinert, April 2020, May 2021.
"""

import time
import calendar
from enum import Enum

import numpy as np

from tools_GESLA import read_GESLA_surge_file, get_GESLA_surge_filename


# Path to the skew surge data file of Reinert et al. (2021)
SKEW_SURGE_FILE = "data/skew_surge_Brest.npy"


class Timeseries(Enum):
    SKEW_SURGE = "skew surge (Reinert et al. 2021)"
    SKEW_SURGE_GESLA = "skew surge (GESLA-2)"


def load_data(city: str, dataset: Timeseries, year_start=None, year_end=None, include_low_tide_surge=False):
    """Get ‘dataset’ for ‘city’ from ‘year_start’ to ‘year_end’.

    The returned object is a dictionary with the following keys:
     - city (string, the same as the input value)
     - name (string, short description of the dataset)
     - full_name (string, long description of the dataset)
     - year_start (int, first valid year of the time series)
     - year_end (int, last valid year of the time series)
     - t (array, time values in s)
     - h (array, height values in cm)

    The returned dataset can be limited to contain only data from
    ‘year_start’ to ‘year_end’.  Each of these two arguments is optional
    and inclusive if given.  They work even if no data for the given
    years exists, in which case the returned values of ‘year_start’ and
    ‘year_end’ will be set to the actual first and last year for which
    there is data within the given period.

    The surge dataset by Reinert et al. (2021) is only available for
    Brest.  It contains both high and low tide skew surge levels, but by
    default, only high tide surge levels are used.  To also include low
    tide surge, set the corresponding parameter to True.  Surge levels
    by Reinert et al. (2021) are relative to the annual mean sea level.

    For the GESLA-2 skew surge dataset by Marcos & Woodworth (2017),
    several cities are available.  This dataset contains only high tide
    surge levels, so include_low_tide_surge must be False.  Surge levels
    in the GESLA-2 dataset contain variations in the mean sea level, in
    particular the mean sea level rise.
    """
    # Initialize the dictionary that is returned
    data = {
        "city": city,
        "name": "",
        "full_name": "",
        "year_start": None,
        "year_end": None,
        "t": None,
        "h": None,
    }
    if dataset == Timeseries.SKEW_SURGE:
        if city != "Brest":
            raise ValueError(
                "only city Brest is available for SKEW_SURGE dataset; "
                "did you want to use SKEW_SURGE_GESLA dataset?"
            )
        data["name"] = "surge (Reinert et al. 2021)"
        data["full_name"] = "skew surge relative to annual mean sea level (Reinert et al. 2021)"
        print("Reading", data["full_name"], "for", data["city"])
        data["t"], data["h"], high_tides = np.load(SKEW_SURGE_FILE)
        high_tides = np.array(high_tides, dtype=bool)
        if include_low_tide_surge:
            print("Using both high and low tide surge levels")
        else:
            data["name"] = "high tide " + data["name"]
            data["full_name"] = "high tide " + data["full_name"]
            data["t"] = data["t"][high_tides]
            data["h"] = data["h"][high_tides]
    elif dataset == Timeseries.SKEW_SURGE_GESLA:
        if include_low_tide_surge:
            raise ValueError("GESLA-2 surge dataset does not contain low tide surge levels.")
        data["name"] = "surge (GESLA-2)"
        data["full_name"] = "skew surge of GESLA-2 (Marcos & Woodworth 2017)"
        filename = get_GESLA_surge_filename(city)
        print("Reading", data["full_name"], "for", data["city"])
        data["t"], data["h"] = read_GESLA_surge_file(filename)
        data["t"] = np.array(data["t"])
        data["h"] = np.array(data["h"])
    else:
        raise ValueError("unknown dataset {} requested".format(dataset))
    # Limit the data to the given range
    if year_start is not None:
        print("Removing data before {}".format(year_start))
        starttime = calendar.timegm((year_start, 1, 1, 0, 0, 0))
        data["h"] = data["h"][data["t"] >= starttime]
        data["t"] = data["t"][data["t"] >= starttime]
    if year_end is not None:
        print("Removing data after {}".format(year_end))
        endtime = calendar.timegm((year_end + 1, 1, 1, 0, 0, 0))
        data["h"] = data["h"][data["t"] < endtime]
        data["t"] = data["t"][data["t"] < endtime]
    # Convert from mm to cm
    data["h"] = data["h"] / 10
    # Get first and last year
    data["year_start"] = time.gmtime(min(data["t"])).tm_year
    data["year_end"] = time.gmtime(max(data["t"])).tm_year
    print("{:9_d} records".format(len(data["t"])))
    return data
