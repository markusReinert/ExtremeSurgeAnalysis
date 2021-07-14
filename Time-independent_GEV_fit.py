"""Make a time-independent GEV fit to monthly maxima (MM).

Written by Markus Reinert, June 2020–July 2021.
"""

import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter

from advanced_GEV_analysis import negative_log_likelihood, XI_THRESHOLD
from advanced_GEV_analysis import get_month_selection, get_year_selection
from tools_errors import format_GEV_parameters, get_error_bounds
from tools_surge import load_data, Timeseries


def GEV(x, mu, sigma, xi):
    """Cumulative distribution function (CDF) of a GEV."""
    if abs(xi) < XI_THRESHOLD:
        # Use a Gumbel distribution
        return np.exp(-np.exp(-((x - mu) / beta)))
    elif xi > 0:
        # For xi > 0, the support of GEV has the lower bound mu-sigma/xi
        x_min = mu - sigma / xi
        y = np.zeros_like(x, dtype=float)
        y[x <= x_min] = 0
        y[x > x_min] = np.exp(-((1 + xi * ((x[x > x_min] - mu) / sigma)) ** (-1 / xi)))
        return y
    else:
        # For xi < 0, the support of GEV has the upper bound mu-sigma/xi
        x_max = mu - sigma / xi
        y = np.zeros_like(x, dtype=float)
        y[x < x_max] = np.exp(-((1 + xi * ((x[x < x_max] - mu) / sigma)) ** (-1 / xi)))
        y[x >= x_max] = 1
        return y


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

# Define a function to compute the return period from a CDF
return_period = lambda P: 1 / (1 - P**(n_months / n_years))

# Fit a GEV to the extreme values
result = optimize.minimize(negative_log_likelihood, [10, 15, -0.1], args=(h_MM,))
if not result.success:
    print("Warning:", result.message)
params = result.x
errors = np.sqrt(np.diag(result.hess_inv))

# Calculate the CDF and the return periods of the fit
x_axis = np.linspace(min(h_MM), max(h_MM), 1000)
P_MM = GEV(x_axis, *params)
T_MM = return_period(P_MM)

# Calculate the 95 % confidence interval of the fit
n_sigma = 1.96
P_MM_bounds = get_error_bounds(GEV, x_axis, params, errors, n_sigma)
T_MM_bounds = [return_period(P_bound) for P_bound in P_MM_bounds]

# Calculate the CDF and the return periods of the empirical data
h_empirical = sorted(h_MM)
P_empirical = (1 + np.arange(n_months)) / (1 + n_months)
T_empirical = return_period(P_empirical)


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
ax.semilogx(T_MM, x_axis, label="GEV fit: " + format_GEV_parameters(params, errors))
assert n_sigma == 1.96, "label 95 % CI is not correct"
ax.fill_betweenx(x_axis, *T_MM_bounds, alpha=0.3, label="95 % confidence interval")

ax.legend()
ax.set_xlim(0.9, 200)
ax.grid(linestyle=":")
ax.xaxis.set_major_formatter(ScalarFormatter())

plt.savefig("results/GEV_fit_{}.png".format(data["city"]))
plt.show()
