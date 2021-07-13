"""Fit a time- and NAO-dependent GEV distribution to monthly maxima (MM).

NAO stands for the index of the North Atlantic Oscillation.
See tools_climate.py for information how to obtain this data.

Written by Markus Reinert, June 2020–July 2021.
"""

import numpy as np
from scipy import optimize, stats

from advanced_GEV_analysis import negative_log_likelihood, Modifiers
from advanced_GEV_analysis import get_month_selection
from advanced_GEV_analysis import compute_amplitude_and_phase, check_significance
from tools_surge import load_data, Timeseries, Subseries
from tools_climate import load_NOAA_data, ClimateIndex


data = load_data("Brest", Timeseries.SKEW_SURGE_GESLA)
data_cli = load_NOAA_data(ClimateIndex.NAO)

# Get the maximum value in every month and its date-time,
# as well as the NAO index of this month
h_MM = []
t_MM = []
cli_MM = []
cli_years = np.array(np.round(data_cli["t_years"], 2), dtype=int)
for year in range(data["year_start"], data["year_end"] + 1):
    for month in range(1, 13):
        sel = get_month_selection(year, month, data["t"])
        if np.any(sel):
            i_max = np.argmax(data["h"][sel])
            value = data["h"][sel][i_max]
            if value < -300:
                # There is an error in the GESLA-2 skew surge
                # dataset for Brest, which we want to catch here.
                print("Neglecting unrealistic outlier of {} m in surge data.".format(value))
                continue
            h_MM.append(value)
            t_MM.append(data["t"][sel][i_max])
            # Check if there is climate data for this year
            if year in cli_years:
                i_year = np.where(cli_years == year)[0]
                # There must be exactly 12 monthly values (with possibly some being NaN)
                assert len(i_year) == 12, "not exactly 12 monthly values of climate index"
                cli_MM.append(data_cli["index"][i_year[month - 1]])
            else:
                cli_MM.append(np.nan)
h_MM = np.array(h_MM)
t_MM = np.array(t_MM)
cli_MM = np.array(cli_MM)
# Replace all NaNs in the climate data by 0s, because adding zero does not change anything
cli_MM[np.isnan(cli_MM)] = 0

# Fit time-independent GEV to the extreme values for a first estimate of the fit parameters
print("")
result = optimize.minimize(negative_log_likelihood, [10, 15, -0.1], args=(h_MM,))
if not result.success:
    print("Warning:", result.message)
params_const = result.x
errors_const = np.sqrt(np.diag(result.hess_inv))
# Compute the (positive) log-likelihood of the model
LL_const = -negative_log_likelihood(params_const, h_MM)
print("\nEstimated GEV parameters for time-independent model:")
for i, name in enumerate(["mu", "sigma", "xi"]):
    print("{:5} = {:7.4f} ± {:.4f}".format(name, params_const[i], errors_const[i]))
print("Parameter values of mu and sigma in cm, xi dimensionless.")


# Fit GEV model with linear trends and annual cycles in the parameters
print("")
args = (
    h_MM,
    t_MM,
    [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION],
    [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION],
)
mu, sigma, xi = params_const
params_initial = [round(mu), 0.1, 1, 1, round(sigma), 0.1, 1, 1, xi]
result = optimize.minimize(negative_log_likelihood, params_initial, args, options={"gtol": 5e-4})
if not result.success:
    print("Warning:", result.message)
params_trends = result.x
errors_trends = np.sqrt(np.diag(result.hess_inv))
LL_trends = -negative_log_likelihood(params_trends, *args)
print("\nEstimated GEV parameters for model with linear trends and annual cycles:")
for i, name in enumerate([
        "mu", "mu_trend", "mu_cos", "mu_sin",
        "sigma", "sigma_trend", "sigma_cos", "sigma_sin", "xi",
]):
    print("{:11} = {:7.4f} ± {:.4f}".format(name, params_trends[i], errors_trends[i]))
print("Trends in cm per century, xi dimensionless, other parameter values in cm.")

print("\nCompare time-dependent GEV model with time-independent model.")
D = 2 * (LL_trends - LL_const)
k = len(params_trends) - len(params_const)
check_significance(D, k)


# Fit GEV model with linear trends, annual cycles, and linear dependencies on climate index
print("")
args = (
    h_MM,
    t_MM,
    [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION, Modifiers.CLIMATE_INDEX],
    [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION, Modifiers.CLIMATE_INDEX],
    None,  # value of shape parameter xi, None means “best fit”, 0.0 means Gumbel distribution
    cli_MM,
)
params_initial = [round(mu), 0.1, 1, 1, 0.1, round(sigma), 0.1, 1, 1, 0.1, xi]
result = optimize.minimize(negative_log_likelihood, params_initial, args, options={"gtol": 5e-4})
if not result.success:
    print("Warning:", result.message)
params_cli = result.x
errors_cli = np.sqrt(np.diag(result.hess_inv))
LL_cli = -negative_log_likelihood(params_cli, *args)
print("\nEstimated GEV parameters for model with linear trends, annual cycles, "
      "and linear {}-dependencies:".format(data_cli["name"]))
for i, name in enumerate([
        "mu", "mu_trend", "mu_cos", "mu_sin", "mu_" + data_cli["name"],
        "sigma", "sigma_trend", "sigma_cos", "sigma_sin", "sigma_" + data_cli["name"], "xi",
]):
    print("{:11} = {:7.4f} ± {:.4f}".format(name, params_cli[i], errors_cli[i]))
print("Trends in cm per century, {} coefficients in cm per unit, "
      "xi dimensionless, other parameter values in cm.".format(data_cli["name"]))

print("\nCompare time-dependent GEV model with linear {}-dependencies "
      "to model without.".format(data_cli["name"]))
D = 2 * (LL_cli - LL_trends)
k = len(params_cli) - len(params_trends)
check_significance(D, k)

# Convert coefficients of cos and sin to amplitude and phase
print("")
for name, i_cos in zip(["mu", "sigma"], [2, 6]):
    print("Annual cycle in {}:".format(name))
    amp, phase, std_amp, std_phase = compute_amplitude_and_phase(
        params_cli[i_cos], params_cli[i_cos + 1],
        errors_cli[i_cos], errors_cli[i_cos + 1],
    )
    print("  Amplitude = ({:.1f} ± {:.1f}) cm".format(amp, std_amp))
    print("  Phase = ({:.0f} ± {:.0f})°".format(np.rad2deg(phase), np.rad2deg(std_phase)))


print(
    "\nNote that surge levels in the GESLA-2 dataset are NOT corrected "
    "for the mean sea level rise, so a trend in mu, comparable to the "
    "mean sea level rise in {}, is expected, when this dataset is used.".format(data["city"])
)
