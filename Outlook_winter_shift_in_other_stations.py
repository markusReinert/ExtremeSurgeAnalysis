"""Analyse extreme surges in different stations with monthly GEV models.

This code implements the analysis presented by Reinert et al. (2021).

Written by Markus Reinert, August 2020 to October 2021.
"""

import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt

from advanced_GEV_analysis import negative_log_likelihood, Modifiers, get_month_selection
from tools_surge import load_data, Timeseries


# Select months of interest (here: March and October)
selected_months = [3, 10]

# Select stations of interest
cities = [
    "Newlyn",
    "Cherbourg",
    "Le Havre",
    "Brest (our data)",
    "Brest (GESLA-2)",
    "La Rochelle",
    "Saint-Jean-de-Luz",
    "Santander",
    "La Coruna",
    "Vigo",
]

# Select the period of interest
year_start = 1950
year_end = 2000


# Define ranges for initial parameter values
# Explanation:
# To fit the time-dependent GEV model to monthly maxima using
# scipy.optimize.minimize, we need to specify an initial value for each
# parameter of the model.  These initial values should be close to
# parameter values of the best fit.  Since the best fit varies from
# station to station and from month to month, there is probably no set
# of initial values that works well for every station.  Therefore, we
# search at first for the best fit over a discrete space of common
# parameter values, and we then use these values as initial values to
# compute the actual best fit.  With the values given here, determining
# the initial values takes between 1 and 2 seconds on a normal PC, which
# we find an acceptable computation time.
mu0_range = np.linspace(5, 40, 5)
mu1_range = np.linspace(-50, 50, 10)
si0_range = np.linspace(5, 25, 5)
si1_range = np.linspace(-50, 50, 10)
xi_range = np.linspace(-0.6, 0.2, 11)
n_combinations = (
    mu0_range.size * mu1_range.size * si0_range.size * si1_range.size * xi_range.size
)
print("Number of initial value combinations:", n_combinations)
t_calculation = 50e-6  # time typically required to compute negative_log_likelihood
print(f"Testing all combinations will take about {n_combinations*t_calculation:.1f} seconds")

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
            options={"gtol": 1e-4},
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
        # Find good initial parameters for the fit
        least_nLL = np.infty
        best_mu0 = fit_params[0]
        best_mu1 = 1.0
        best_si0 = fit_params[4]
        best_si1 = 1.0
        best_xi = fit_params[-1]
        for mu0 in mu0_range:
            for mu1 in mu1_range:
                for si0 in si0_range:
                    for si1 in si1_range:
                        for xi in xi_range:
                            nLL = negative_log_likelihood(
                                [mu0, mu1, si0, si1, xi],
                                h_MM[month],
                                t_MM[month],
                                [Modifiers.LINEAR_TREND],
                                [Modifiers.LINEAR_TREND],
                            )
                            if nLL < least_nLL:
                                least_nLL = nLL
                                best_mu0 = mu0
                                best_mu1 = mu1
                                best_si0 = si0
                                best_si1 = si1
                                best_xi = xi
        # Fit GEV model with linear trends in the parameters for this month
        with np.errstate(invalid="ignore"):
            # Do not warn when the function value infinity occurs in the optimisation
            result_month = optimize.minimize(
                negative_log_likelihood,
                [best_mu0, best_mu1, best_si0, best_si1, best_xi],
                args=(
                    h_MM[month], t_MM[month], [Modifiers.LINEAR_TREND], [Modifiers.LINEAR_TREND]
                ),
                options={"gtol": 1e-4},
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
        cities,
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
