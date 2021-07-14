"""Analyse extreme surges month by month with time-dependent GEV distributions.

This code implements the “monthly analysis” (method 2) of Reinert et
al. (2021).

Written by Markus Reinert, July 2020–July 2021.
"""

from calendar import timegm

import numpy as np
from scipy import optimize, special
from matplotlib import pyplot as plt

from advanced_GEV_analysis import negative_log_likelihood, Modifiers, get_month_selection
from advanced_GEV_analysis import time_dependent_GEV_parameters, TREND_SCALE
from tools_surge import load_data, Timeseries


months = np.arange(1, 13)
month_names = [
    "Jan", "Feb", "Mar", "April", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"
]
month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def GEV_mean(mu, sigma, xi):
    """Calculate the mean of a GEV distribution.

    The arguments mu, sigma, and xi are the location, scale, and shape
    parameter of the GEV distribution, respectively.
    """
    if xi >= 1:
        return np.inf
    elif xi == 0:
        return mu + sigma * np.euler_gamma
    else:
        g1 = special.gamma(1 - xi)
        return mu + sigma * (g1 - 1) / xi


def start_year_in_july(data):
    month = 7
    return [*data[month-1:], *data[:month-1]]


def select_winter(monthly_data):
    return start_year_in_july(monthly_data)[3:-3]


transp_area = 0.3
n_sigma = 1.96


# Choose a period and load the surge data
data = load_data("Brest", Timeseries.SKEW_SURGE_GESLA, year_start=1950, year_end=2000)


# Calculate monthly GEV fits with linear trends
parameter_names = ["mu", "mu_trend", "sigma", "sigma_trend", "xi"]
fit_params_monthly = []
fit_errors_monthly = []
for month in months:
    print("\nMonth:", month)
    # Get the maximum value for this month in every year
    t_MM = []
    h_MM = []
    for year in range(data["year_start"], data["year_end"] + 1):
        # Uncomment this to check the influence of the Great Storm of 1987
        # if year == 1987:
        #     print("Skipping year 1987")
        #     continue
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
    t_MM = np.array(t_MM)
    h_MM = np.array(h_MM)
    print("number of years =", len(t_MM))

    # Fit time-independent GEV for a first estimate of the fit parameters
    result_no_trend = optimize.minimize(
        negative_log_likelihood, [20, 15, -0.1], args=(h_MM,), options={"gtol": 5e-5}
    )
    if not result_no_trend.success:
        print("Warning:", result_no_trend.message, "(fit without trends)")
    mu, sigma, xi = result_no_trend.x

    # Fit GEV model with linear trends in the parameters
    result_both_trends = optimize.minimize(
        negative_log_likelihood,
        [mu, 1, sigma, 1, xi],
        args=(h_MM, t_MM, [Modifiers.LINEAR_TREND], [Modifiers.LINEAR_TREND]),
        options={"gtol": 5e-5},
    )
    if not result_both_trends.success:
        print("Warning:", result_both_trends.message, "(fit with both trends)")

    # Print and save parameter estimates and their standard errors
    fit_params = result_both_trends.x
    fit_errors = np.sqrt(np.diag(result_both_trends.hess_inv))
    for i, name in enumerate(parameter_names):
        print("{:11} = {:7.4f} ± {:.4f}".format(name, fit_params[i], fit_errors[i]))
    fit_params_monthly.append(fit_params)
    fit_errors_monthly.append(fit_errors)


# Show monthly variation of all fit parameters
fig, axs = plt.subplots(nrows=len(parameter_names), sharex=True)
fig.suptitle("{}: {} to {}".format(data["city"], data["year_start"], data["year_end"]))
for i, ax in enumerate(axs):
    ax.set_title(parameter_names[i])
    ax.plot(month_names, [params[i] for params in fit_params_monthly], "o-")
    ax.fill_between(
        month_names,
        [
            params[i] + n_sigma * errors[i]
            for params, errors in zip(fit_params_monthly, fit_errors_monthly)
        ],
        [
            params[i] - n_sigma * errors[i]
            for params, errors in zip(fit_params_monthly, fit_errors_monthly)
        ],
        alpha=transp_area,
    )
    ax.grid()
plt.show(block=False)


fig, ax = plt.subplots()

ax.set_title("Mean extreme surge per month in {}".format(data["city"]))
ax.set_ylabel("Surge level in cm")

# To distribute the fit uncertainty equally between the start year and
# the end year, shift the Epoch in the middle between start and end
t_start = timegm((data["year_start"], 1, 1, 0, 0, 0))
t_end = timegm((data["year_end"], 12, 1, 0, 0, 0))
t_offset = t_start + (t_end - t_start) / 2

for year in [data["year_start"], data["year_end"]]:
    # Compute the mean of the GEV distribution of every month and its standard error
    monthly_means = np.zeros(len(months))
    monthly_means_std = np.zeros(len(months))
    for i, month in enumerate(months):
        time_secs = timegm((year, month, 1, 0, 0, 0))
        mu, sigma, xi = time_dependent_GEV_parameters(
            fit_params_monthly[i],
            time_secs,
            [Modifiers.LINEAR_TREND],
            [Modifiers.LINEAR_TREND],
        )
        monthly_means[i] = GEV_mean(mu, sigma, xi)
        monthly_means_std[i] = np.sqrt(
            fit_errors_monthly[i][0]**2
            + (fit_errors_monthly[i][1] * (time_secs - t_offset) / TREND_SCALE)**2
            + (
                fit_errors_monthly[i][2]**2
                + (fit_errors_monthly[i][3] * (time_secs - t_offset) / TREND_SCALE)**2
            ) * ((special.gamma(1 - xi) - 1) / xi)**2
            + (
                fit_errors_monthly[i][4] * sigma / xi**2 * (
                    1 - special.gamma(1-xi) * (1 + xi * special.polygamma(0, 1-xi))
                )
            )**2
        )
    ax.plot(
        start_year_in_july(month_names),
        start_year_in_july(monthly_means),
        "o-",
        label=str(year),
    )
    ax.fill_between(
        start_year_in_july(month_names),
        start_year_in_july(monthly_means + n_sigma * monthly_means_std),
        start_year_in_july(monthly_means - n_sigma * monthly_means_std),
        alpha=transp_area,
    )

ax.legend()
ax.grid(linestyle=":")

plt.show()


# Calculate the date of mid-winter from the mu-values
print("")
print("Analysis of winter:")
winter_months = select_winter(months)
winter_month_days = select_winter(month_days)
winter_month_names = select_winter(month_names)
for year in [data["year_start"], data["year_end"]]:
    # Calculate mu at the beginning of each winter month
    winter_month_mus = [
        time_dependent_GEV_parameters(
            fit_params_monthly[month-1],
            timegm((year, month, 1, 0, 0, 0)),
            [Modifiers.LINEAR_TREND],
            [Modifiers.LINEAR_TREND],
        )[0]
        for month in winter_months
    ]
    print("*", year)
    print(" ", "winter mu values:")
    print(" ", "  ".join(winter_month_names))
    print(" ", " ".join(["{:4.1f}".format(mu) for mu in winter_month_mus]))
    winter_days = 1 + np.arange(sum(winter_month_days))
    weights = []
    for days, mu_value in zip(winter_month_days, winter_month_mus):
        weights.extend([mu_value] * days)
    mid_winter_decimal = np.average(winter_days, weights=weights)
    i_month = 0
    while mid_winter_decimal > winter_month_days[i_month]:
        mid_winter_decimal -= winter_month_days[i_month]
        i_month += 1
    print(" ", "weighted mid-winter: {day:2.0f} {month}".format(
        day=mid_winter_decimal, month=winter_month_names[i_month]
    ))
