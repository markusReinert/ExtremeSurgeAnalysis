"""Make a time-independent GEV fit to monthly maxima (MM).

Written by Markus Reinert, June 2020–March 2022.
"""

import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter

from advanced_GEV_analysis import negative_log_likelihood, get_year_selection, get_month_selection
from advanced_GEV_analysis import GEV_return_level, GEV_standard_error, format_GEV_parameters
from tools_surge import load_data, Timeseries


data = load_data("Brest", Timeseries.SKEW_SURGE_GESLA)

# Get the maximum value in every month
h_MM = []
n_years = 0
for year in range(data["year_start"], data["year_end"] + 1):
    sel_year = get_year_selection(year, data["t"])
    if any(sel_year):
        n_years += 1
        for month in range(1, 13):
            sel = get_month_selection(year, month, data["t"])
            if np.any(sel):
                value = max(data["h"][sel])
                if value < -300:
                    # There is an error in the GESLA-2 skew surge
                    # dataset for Brest, which we want to catch here.
                    print("Neglecting unrealistic outlier of {} m.".format(value))
                else:
                    h_MM.append(value)
h_MM = np.array(h_MM)
n_months = len(h_MM)

# Fit a GEV to the extreme values
result = optimize.minimize(negative_log_likelihood, [10, 15, -0.1], args=(h_MM,))
if not result.success:
    print("Warning:", result.message)
params = result.x
covars = result.hess_inv
errors = np.sqrt(np.diag(covars))

# Calculate the graph and the standard error of the fitted GEV model
t_axis = np.logspace(0, 3, 10_000)[1:]  # exclude t = 1 to avoid value -infinity
h_model = GEV_return_level(t_axis, *params, values_per_year=n_months/n_years)
h_std_error = GEV_standard_error(t_axis, *params, covars, values_per_year=n_months/n_years)

# Calculate the CDF and the return periods of the empirical data
h_empirical = sorted(h_MM)
P_empirical = (1 + np.arange(n_months)) / (1 + n_months)
T_empirical = 1 / (1 - P_empirical**(n_months / n_years))


fig, ax = plt.subplots()

ax.set_title(
    "GEV fit to monthly surge maxima in {} from {} to {}".format(
        data["city"], data["year_start"], data["year_end"]
    ),
    weight="bold",
)
ax.set_xlabel("Return period in years")
ax.set_ylabel("Return level in cm")

ax.semilogx(
    T_empirical, h_empirical, "k.",
    label="Empirical return periods ({} data points)".format(n_months),
)
ax.semilogx(t_axis, h_model, label="GEV fit: " + format_GEV_parameters(params, errors))
ax.fill_between(
    t_axis,
    h_model + 1.96 * h_std_error,
    h_model - 1.96 * h_std_error,
    alpha=0.3,
    label="95 % confidence interval",
)

ax.legend()
ax.set_xlim(0.9, 200)
ax.grid(linestyle=":")
ax.xaxis.set_major_formatter(ScalarFormatter())

plt.savefig("results/GEV_fit_{}.png".format(data["city"]))
plt.show()
