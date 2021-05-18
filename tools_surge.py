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


class Subseries(Enum):
    HIGH_TIDE_SURGE = "high tide surge"
    LOW_TIDE_SURGE = "low tide surge"


def load_data(city: str, dataset: Timeseries, year_start=None, year_end=None, qualifiers=[]):
    """Get ‘dataset’ for ‘city’ from ‘year_start’ to ‘year_end’ after applying ‘qualifiers’.

    The returned object is a dictionary with the following keys:
     - city (string, the same as the input value)
     - name (string, short description of the dataset)
     - full_name (string, long description of the dataset)
     - year_start (int, first valid year of the timeseries)
     - year_end (int, last valid year of the timeseries)
     - t (array, time values in s)
     - h (array, height values in cm)

    The returned dataset can be limited to contain only data from
    ‘year_start’ to ‘year_end’.  Each of these two arguments is optional
    and inclusive if given.  They work even if no data for the given
    years exists, in which case the returned values of ‘year_start’ and
    ‘year_end’ will be set to the actual first and last year for which
    there is data within the given period

    Note that SKEW_SURGE_GESLA contains only high tide surge and
    includes variations in the mean sea level.  On the other hand,
    SKEW_SURGE is only available for Brest, uses a yearly tidal
    prediction so that the mean sea level rise is not part of the surge,
    and can contain high tide surge or low tide surge or both.  Use the
    argument ‘qualifiers’ to choose between high or low tide surge.
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
        data["full_name"] = "skew surge from yearly predictions (Reinert et al. 2021)"
        print("Reading", data["full_name"], "for", data["city"])
        data["t"], data["h"], high_tides = np.load(SKEW_SURGE_FILE)
        high_tides = np.array(high_tides, dtype=bool)
        if Subseries.HIGH_TIDE_SURGE in qualifiers and Subseries.LOW_TIDE_SURGE in qualifiers:
            # High and low tide surge together give all surge levels
            pass
        elif Subseries.HIGH_TIDE_SURGE in qualifiers:
            print("Restricting to high tide surge")
            data["name"] = "high tide " + data["name"]
            data["full_name"] = "high tide " + data["full_name"]
            data["t"] = data["t"][high_tides]
            data["h"] = data["h"][high_tides]
        elif Subseries.LOW_TIDE_SURGE in qualifiers:
            print("Restricting to low tide surge")
            data["name"] = "low tide " + data["name"]
            data["full_name"] = "low tide " + data["full_name"]
            data["t"] = data["t"][~high_tides]
            data["h"] = data["h"][~high_tides]
        else:
            # No restriction to high or low tide, that means all surge levels
            pass
    elif dataset == Timeseries.SKEW_SURGE_GESLA:
        if qualifiers:
            raise ValueError("GESLA-2 surge dataset does not take any qualifiers.")
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
