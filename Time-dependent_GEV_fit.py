"""Fit a GEV distribution with time-dependent parameters to monthly maxima (MM).

Written by Markus Reinert, June 2020–July 2021.
"""

import numpy as np
from scipy import optimize

from advanced_GEV_analysis import negative_log_likelihood, Modifiers
from advanced_GEV_analysis import get_month_selection
from advanced_GEV_analysis import compute_amplitude_and_phase, check_significance
from tools_surge import load_data, Timeseries, Subseries


data = load_data("Brest", Timeseries.SKEW_SURGE_GESLA)

# Get the maximum value in every month and its date-time
h_MM = []
t_MM = []
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
h_MM = np.array(h_MM)
t_MM = np.array(t_MM)


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


# Fit GEV model with annual cycles in the parameters
print("")
args = (h_MM, t_MM, [Modifiers.SEASONAL_OSCILLATION], [Modifiers.SEASONAL_OSCILLATION])
mu, sigma, xi = params_const
params_initial = [round(mu), 1, 1, round(sigma), 1, 1, xi]
result = optimize.minimize(negative_log_likelihood, params_initial, args)
if not result.success:
    print("Warning:", result.message)
params_annual = result.x
errors_annual = np.sqrt(np.diag(result.hess_inv))
LL_annual = -negative_log_likelihood(params_annual, *args)
print("\nEstimated GEV parameters for model with annual cycles:")
for i, name in enumerate(["mu", "mu_cos", "mu_sin", "sigma", "sigma_cos", "sigma_sin", "xi"]):
    print("{:9} = {:7.4f} ± {:.4f}".format(name, params_annual[i], errors_annual[i]))
print("Parameter values in cm, except for xi (dimensionless).")

print("\nCompare GEV model with annual cycles to time-independent GEV model.")
D = 2 * (LL_annual - LL_const)
k = len(params_annual) - len(params_const)
check_significance(D, k)


# Fit GEV model with annual cycles and linear trends in the parameters
print("")
args = (
    h_MM,
    t_MM,
    [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION],
    [Modifiers.LINEAR_TREND, Modifiers.SEASONAL_OSCILLATION],
)
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

print("\nCompare GEV model with linear trends and annual cycles "
      "to GEV model with only annual cycles.")
D = 2 * (LL_trends - LL_annual)
k = len(params_trends) - len(params_annual)
check_significance(D, k)
print(
    "Note that surge levels in the GESLA-2 dataset are NOT corrected "
    "for the mean sea level rise, so a trend in mu, comparable to the "
    "mean sea level rise in {}, is expected, when this dataset is used.".format(data["city"])
)

# Convert coefficients of cos and sin to amplitude and phase
print("")
for name, i_cos in zip(["mu", "sigma"], [2, 6]):
    print("Annual cycle in {}:".format(name))
    amp, phase, std_amp, std_phase = compute_amplitude_and_phase(
        params_trends[i_cos], params_trends[i_cos + 1],
        errors_trends[i_cos], errors_trends[i_cos + 1],
    )
    print("  Amplitude = ({:.1f} ± {:.1f}) cm".format(amp, std_amp))
    print("  Phase = ({:.0f} ± {:.0f})°".format(np.rad2deg(phase), np.rad2deg(std_phase)))
