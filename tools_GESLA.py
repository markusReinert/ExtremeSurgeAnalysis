"""Handling of GESLA-2 skew surge data files.

Surge data files can be found on https://www.gesla.org/
at the end of the reference “Marcos, M. and Woodworth, P.L. 2017”
(last checked on 18 May 2021).

Extract those that you are interested in and place them in the folder
specified as ‘GESLA_DIR’.  Make sure that their names and the
corresponding city are correctly implemented in the
‘get_GESLA_surge_filename’ method.

Written by Markus Reinert, May 2020/2021.
"""

import os
import calendar


GESLA_DIR = "data/skew_surge_GESLA/"


def get_GESLA_surge_filename(city):
    """Get the path to the file containing GESLA-2 surge data for ‘city’."""
    if city == "Brest":
        filename = "brest-brest-france-refmar_SkewSurges.txt"
    elif city == "La Rochelle":
        filename = "la_rochelle_la_palli-la_rochelle_la_palli-france-refmar_SkewSurges.txt"
    elif city == "Saint-Jean-de-Luz":
        filename = "saint_jean_de_luz_so-saint_jean_de_luz_so-france-refmar_SkewSurges.txt"
    elif city == "Cherbourg":
        filename = "cherbourg-cherbourg-france-refmar_SkewSurges.txt"
    elif city == "Le Havre":
        filename = "le_havre-le_havre-france-refmar_SkewSurges.txt"
    elif city == "Newlyn":
        # Choose one of the following files
        # filename = "newlyn,_cornwall-294a-united_kingdom-uhslc_SkewSurges.txt"  # 1915-2010, 2 MB
        # filename = "newlyn-newlyn-glossdm-bodc_SkewSurges.txt"  # 1916-1944, 600 kB
        filename = "newlyn-p001-uk-bodc_SkewSurges.txt"  # 1915-2014, 1.7 MB
    elif city == "Vigo":
        filename = "vigo-vigo-spain-ieo_SkewSurges.txt"
    elif city == "La Coruna":
        filename = "la_coruna-830a-spain-uhslc_SkewSurges.txt"
    elif city == "Santander":
        filename = "santander-sant-spain-ieo_SkewSurges.txt"
    else:
        raise ValueError("unknown city: " + repr(city))
    return os.path.join(GESLA_DIR, filename)


def read_GESLA_surge_file(filename):
    """Read a file with skew surge data from the GESLA-2 dataset.

    Return two lists of integers, containing the
    time in seconds and the height in mm, respectively.

    The input file must contain lines in the format
    ‘year month day hour surge’, where the surge is given in meter.
    Lines starting with "#" are skipped.
    """
    time_values = []
    surge_values = []
    with open(filename) as f:
        for i, line in enumerate(f.readlines()):
            if line.startswith("#"):
                continue
            year, month, day, hour, surge = line.split()
            # Convert the date-time to seconds since the epoch (1970-01-01)
            time_values.append(
                int(calendar.timegm((int(year), int(month), int(day), int(hour), 0, 0)))
            )
            surge_values.append(int(round(1000 * float(surge))))
    return time_values, surge_values
