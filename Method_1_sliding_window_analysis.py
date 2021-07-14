"""Fit a time-dependent GEV distribution to monthly maxima in a sliding window.

This code implements the “sliding window analysis” (method 1) of Reinert
et al. (2021).

Written by Markus Reinert, October 2020–July 2021.
"""

from datetime import date
from calendar import timegm

import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt

from advanced_GEV_analysis import negative_log_likelihood, Modifiers, get_month_selection
from advanced_GEV_analysis import compute_amplitude_and_phase
from tools_surge import load_data, Timeseries, Subseries


data = load_data("Brest", Timeseries.SKEW_SURGE_GESLA)

# Get time (t), height (h), and year (y) of the monthly maxima (MM)
t_MM = []
h_MM = []
y_MM = []
for year in range(data["year_start"], data["year_end"] + 1):
    for month in range(1, 13):
        sel = get_month_selection(year, month, data["t"])
        if np.any(sel):
            i_max = np.argmax(data["h"][sel])
            value = data["h"][sel][i_max]
            if value < -300:
                # There is an error in the GESLA-2 skew surge
                # dataset for Brest, which we want to catch here.
                print("Neglecting unrealistic outlier of {} m.".format(value))
            else:
                h_MM.append(value)
                t_MM.append(data["t"][sel][i_max])
                y_MM.append(year)
t_MM = np.array(t_MM)
h_MM = np.array(h_MM)
y_MM = np.array(y_MM)


# Choose the size of the sliding window (in years)
window_years = 30
# Choose the minimum number of months required in each window
minimum_months = 6 * window_years

# Compute start point, center point, and end point of the sliding windows
start_years = np.arange(data["year_start"] - window_years, data["year_end"] + 1)
center_years = start_years + window_years / 2 - 1
end_years = start_years + window_years - 1

# Choose initial parameters for the fit in the following order:
# mu, mu_trend, mu_cos, mu_sin, sigma, sigma_trend, sigma_cos, sigma_sin, xi
initial_parameters = [18, 0.1, 1, 1, 15, 0.1, 1, 1, -0.1]
fit_params = np.zeros((len(initial_parameters), len(start_years)))
fit_errors = np.zeros_like(fit_params)
for i_year, start_year in enumerate(start_years):
    end_year = end_years[i_year]
    print("Sliding window from {} to {}".format(start_year, end_year))
    sel_MM = (y_MM >= start_year) & (y_MM <= end_year)
    available_months = np.count_nonzero(sel_MM)
    if available_months >= minimum_months:
        # Do not warn when the function value infinity occurs in the optimisation
        with np.errstate(invalid="ignore"):
            result = optimize.minimize(
                negative_log_likelihood,
                initial_parameters,
                args=(
                    h_MM[sel_MM],
                    t_MM[sel_MM],
                    [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION],
                    [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION],
                ),
                options={"gtol": 1e-4},
            )
        if not result.success:
            print("Warning:", result.message)
        fit_params[:, i_year] = result.x
        fit_errors[:, i_year] = np.sqrt(np.diag(result.hess_inv))
    else:
        print("not enough months:", available_months, "<", minimum_months)
        fit_params[:, i_year] = np.nan
        fit_errors[:, i_year] = np.nan

mu_amp, mu_phase, mu_std_amp, mu_std_phase = compute_amplitude_and_phase(
    fit_params[2], fit_params[3], fit_errors[2], fit_errors[3]
)


fig, ax = plt.subplots()
ax.set_title("Timing of the extreme surge season in {}".format(data["city"]))

ax.plot(center_years, mu_phase, label="method 1")
ax.fill_between(
    center_years,
    mu_phase + 1.96 * mu_std_phase,
    mu_phase - 1.96 * mu_std_phase,
    alpha=0.3,
    label="95 % confidence interval",
)
ax.legend()

length_year = 24 * 3600 * 365.2425
yticks = np.array([
    timegm((1969, 12, 1, 0, 0, 0)),
    timegm((1970,  1, 1, 0, 0, 0)),
    timegm((1970,  2, 1, 0, 0, 0)),
    timegm((1970,  3, 1, 0, 0, 0)),
])
ax.yaxis.set_ticks(2 * np.pi * yticks / length_year)
ax.yaxis.set_ticklabels([date.fromtimestamp(t).strftime("%d %b") for t in yticks])
ax.set_ylim(2 * np.pi * min(yticks) / length_year, 2 * np.pi * max(yticks) / length_year)
ax.set_xlim(data["year_start"], data["year_end"])
ax.grid(linestyle=":")

plt.show()
