"""Analyse extreme surges in different stations with monthly GEV models.

This code implements the analysis presented in the Discussion of the
manuscript by Reinert et al. (2021).

Written by Markus Reinert, August 2020–July 2021.
"""

import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt

from advanced_GEV_analysis import negative_log_likelihood, Modifiers, get_month_selection
from tools_surge import load_data, Timeseries


# Select months of interest (here: March and October)
selected_months = [3, 10]

# Select stations of interest and specify the best possible error tolerance
# for the fit over the full year and over each selected month
cities = {
    "Newlyn": {"year": 5e-5, 3: 1e-5, 10: 5e-5},
    "Cherbourg": {"year": 5e-5, 3: 1e-5, 10: 1e-5},
    "Le Havre": {"year": 5e-4, 3: 4e-5, 10: 1e-4},
    "Brest (our data)": {"year": 1e-5, 3: 5e-5, 10: 5e-5},
    "Brest (GESLA-2)": {"year": 5e-5, 3: 1e-4, 10: 1e-5},
    "La Rochelle": {"year": 1e-5, 3: 1e-5, 10: 1e-5},
    "Saint-Jean-de-Luz": {"year": 5e-5, 3: 5e-5, 10: 5e-5},
    "Santander": {"year": 5e-5, 3: 1e-5, 10: 1e-5},
    "La Coruna": {"year": 5e-5, 3: 1e-5, 10: 1e-5},
    "Vigo": {"year": 1e-4, 3: 1e-5, 10: 5e-5},
}

# Select the period of interest
year_start = 1950
year_end = 2000


# Compute for each station the linear trend in each month of interest
# and over the whole year, as well as their standard errors
linear_trends = {month: [] for month in selected_months}
linear_trends["year"] = []
linear_trends_std = {month: [] for month in selected_months}
linear_trends_std["year"] = []
for city in cities:
    print("")
    # Load the surge dataset of this city for the period of interest
    if city == "Brest (our data)":
        data = load_data("Brest", Timeseries.SKEW_SURGE, year_start, year_end)
    elif city == "Brest (GESLA-2)":
        data = load_data("Brest", Timeseries.SKEW_SURGE_GESLA, year_start, year_end)
    else:
        data = load_data(city, Timeseries.SKEW_SURGE_GESLA, year_start, year_end)
    # Get (date-)time and height of all monthly maxima,
    # and sort those of interest by month
    t_MM = {month: [] for month in selected_months}
    t_MM["year"] = []
    h_MM = {month: [] for month in selected_months}
    h_MM["year"] = []
    for year in range(data["year_start"], data["year_end"] + 1):
        for month in range(1, 13):
            sel = get_month_selection(year, month, data["t"])
            if np.any(sel):
                i_max = np.argmax(data["h"][sel])
                t_MM["year"].append(data["t"][sel][i_max])
                h_MM["year"].append(data["h"][sel][i_max])
                if month in selected_months:
                    t_MM[month].append(data["t"][sel][i_max])
                    h_MM[month].append(data["h"][sel][i_max])
    for k in t_MM:
        t_MM[k] = np.array(t_MM[k])
        h_MM[k] = np.array(h_MM[k])
    # Fit a GEV model with annual cycles and linear trends
    # to all monthly maxima to extract the overall trend
    with np.errstate(invalid="ignore"):
        # Do not warn when the function value infinity occurs in the optimisation
        result_full_year = optimize.minimize(
            negative_log_likelihood,
            [20, 1, 1, 1, 15, 1, 1, 1, -0.1],
            args=(
                h_MM["year"],
                t_MM["year"],
                [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION],
                [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION],
            ),
            options={"gtol": cities[city]["year"]},
        )
    if not result_full_year.success:
        print("Warning:", result_full_year.message)
    fit_params = result_full_year.x
    fit_errors = np.sqrt(np.diag(result_full_year.hess_inv))
    for i, name in enumerate([
            "mu", "mu_trend", "mu_cos", "mu_sin", "sigma", "si_trend", "si_cos", "si_sin", "xi"
    ]):
        print("{:8} = {:7.4f} ± {:.4f}".format(name, fit_params[i], fit_errors[i]))
    linear_trends["year"].append(fit_params[1])
    linear_trends_std["year"].append(fit_errors[1])
    # Compute the linear trends in the months of interest
    for month in selected_months:
        print("Month:", month)
        print("#years =", t_MM[month].size)
        # Fit time-independent GEV model for this month to obtain good initial parameters
        with np.errstate(invalid="ignore"):
            # Do not warn when the function value infinity occurs in the optimisation
            result_timeindep = optimize.minimize(
                negative_log_likelihood,
                [fit_params[0], fit_params[4], fit_params[-1]],
                args=(h_MM[month],),
                options={"gtol": cities[city][month]},
            )
        if not result_timeindep.success:
            print("Warning:", result_timeindep.message, "(time-independent model)")
        # Fit GEV model with linear trends in the parameters for this month
        with np.errstate(invalid="ignore"):
            # Do not warn when the function value infinity occurs in the optimisation
            result_month = optimize.minimize(
                negative_log_likelihood,
                [result_timeindep.x[0], 1, result_timeindep.x[1], 1, result_timeindep.x[2]],
                args=(
                    h_MM[month], t_MM[month], [Modifiers.LINEAR_TREND], [Modifiers.LINEAR_TREND]
                ),
                options={"gtol": cities[city][month]},
            )
        if not result_month.success:
            print("Warning:", result_month.message)
        fit_params_month = result_month.x
        fit_errors_month = np.sqrt(np.diag(result_month.hess_inv))
        for i, name in enumerate(["mu", "mu_trend", "sigma", "si_trend", "xi"]):
            print("{:8} = {:7.4f} ± {:.4f}".format(
                name, fit_params_month[i], fit_errors_month[i]
            ))
        linear_trends[month].append(fit_params_month[1])
        linear_trends_std[month].append(fit_errors_month[1])


fig, ax = plt.subplots()
ax.set_title(
    "Relative linear trend of the location parameter for the period {} to {}".format(
        year_start, year_end
    )
)
ax.set_ylabel("Difference of monthly and yearly trends [cm/century]")
for month in selected_months:
    ax.errorbar(
        cities.keys(),
        [
            linear_trends[month][i_city] - linear_trends["year"][i_city]
            for i_city in range(len(cities))
        ],
        [
            np.sqrt(linear_trends_std[month][i_city]**2 + linear_trends_std["year"][i_city]**2)
            for i_city in range(len(cities))
        ],
        fmt="o",
        label="Month: {}".format(month),
    )
ax.legend()
ax.grid(linestyle=":")
plt.show()
