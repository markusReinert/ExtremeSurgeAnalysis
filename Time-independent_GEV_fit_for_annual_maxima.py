"""Make a time-independent GEV fit to annual maxima (AM).

Written by Markus Reinert, June 2020–March 2022.
"""

import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter

from advanced_GEV_analysis import negative_log_likelihood, get_year_selection
from advanced_GEV_analysis import GEV_return_level, GEV_standard_error, format_GEV_parameters
from tools_surge import load_data, Timeseries


data = load_data("Brest", Timeseries.SKEW_SURGE_GESLA)

# Get the maximum value in every calendar year
h_AM = []
for year in range(data["year_start"], data["year_end"] + 1):
    sel = get_year_selection(year, data["t"])
    if any(sel):
        h_AM.append(max(data["h"][sel]))
h_AM = np.array(h_AM)
n_years = len(h_AM)

# Fit a GEV to the extreme values
result = optimize.minimize(negative_log_likelihood, [40, 15, -0.1], args=(h_AM,))
if not result.success:
    print("Warning:", result.message)
params = result.x
covars = result.hess_inv
errors = np.sqrt(np.diag(covars))

# Calculate the graph and the standard error of the fitted GEV model
t_axis = np.logspace(0, 3, 100_000)[1:]  # exclude t = 1 to avoid value -infinity
h_model = GEV_return_level(t_axis, *params)
h_std_error = GEV_standard_error(t_axis, *params, covars)

# Calculate the CDF and the return periods of the empirical data
h_empirical = sorted(h_AM, reverse=True)
P_empirical = (1 + np.arange(n_years)) / (1 + n_years)
T_empirical = 1 / P_empirical


fig, ax = plt.subplots()

ax.set_title(
    "GEV fit to annual surge maxima in {} from {} to {}".format(
        data["city"], data["year_start"], data["year_end"]
    ),
    weight="bold",
)
ax.set_xlabel("Return period in years")
ax.set_ylabel("Return level in cm")

ax.semilogx(
    T_empirical, h_empirical, "k.",
    label="Empirical return periods ({} data points)".format(n_years),
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
# Limit the y-extent to the observed range +/- 5%
delta_h = (h_AM.max() - h_AM.min()) * 0.05
ax.set_ylim(h_AM.min() - delta_h, h_AM.max() + delta_h)
ax.grid(linestyle=":")
ax.xaxis.set_major_formatter(ScalarFormatter())

plt.savefig("results/GEV_fit_{}_annual.png".format(data["city"]))
plt.show()
