"""Methods for statistical analysis with a time-dependent GEV distribution.

The main reference for this code is the book “An Introduction to
Statistical Modeling of Extreme Values” by Stuart Coles (2001).

The heart of this library is the function negative_log_likelihood, which
is intended to be used for a maximum likelihood fit of a (possibly
time-dependent) GEV function to extreme values.  For convenience, the
function time_dependent_GEV_parameters is provided; it takes almost the
same arguments and returns the three GEV parameters mu, sigma, and xi
corresponding to the given points in time.

Example usage for a time-dependent GEV fit with a linear trend and an
annual oscillation in both parameters (location and scale):

import numpy as np
from scipy import optimize

initial_parameters = [
    # mu
    20, 1, 1, 1,                # average, trend, cosine, sine
    # sigma
    10, 1, 1, 1,                # average, trend, cosine, sine
    # xi
    0.1,
]
result = optimize.minimize(
    negative_log_likelihood,
    initial_parameters,
    args=(
        h_array,                # the observed extreme values
        t_array,                # the times of the events
        [AGA.Modifiers.LINEAR_TREND, AGA.Modifiers.SEASONAL_OSCILLATION],
        [AGA.Modifiers.LINEAR_TREND, AGA.Modifiers.SEASONAL_OSCILLATION],
    ),
)
if not result.success:
    print("Warning:", result.message)
# Get the fitted parameters (array of the same size as initial_parameters)
fit_params = result_GEV.x
# Get the standard deviation corresponding to each fitted parameter
fit_errors = np.sqrt(np.diag(result.hess_inv))


For a time-independent GEV fit, it simplifies to:
initial_parameters = [20, 10, 0.1]  # mu, sigma, xi
result = optimize.minimize(
    negative_log_likelihood, initial_parameters, args=(h_array,)
)


According to Coles (2001), the shape parameter should be bigger than
-0.5 for the maximum likelihood fit to work.  If the fit fails, one can
try to impose this as a strict lower bound by giving the argument
bounds=[
    (None, None),
    ...,
    (None, None),
    (-0.5, None),
]
with the appropriate number of “(None, None)”s (number of
initial_parameters minus 1) to the minimize function.  This changes
automatically the fit method, in which case it is necessary to use
fit_errors = np.sqrt(np.diag(result.hess_inv.todense()))
to obtain the standard deviations.

If the fit fails because the desired error was not necessarily achieved
due to precision loss, then look at the Jacobian (result.jac) and impose
an appropriate error tolerance by giving the argument, e.g.,
options={"gtol": 1e-4}
to the minimize function.

Possible extensions and improvements:
 -- using two climate indices at the same time
 -- make t_array optional if the model depends only on a climate index
 -- The performance of a maximum likelihood fit could be improved by
    implementing the negative_log_likelihood function in a different
    way: Instead of passing a t_array and modifiers for mu and sigma,
    pass directly the arrays that are multiplied by the parameters and
    added together to form the parameters.  This can be, currently, the
    linear t-data, the t-data in a cos and in a sin, and the climate
    indices.  To implement this conveniently, write a new function that
    takes as arguments t_array, climate_index_array, mu_modifiers,
    sigma_modifiers; and returns the arrays that you have to give as
    arrays to the new negative_log_likelihood.  This would also make the
    implementation of the climate index smoother.  In a first approach,
    one can limit this to the case of equal modifiers for mu and sigma,
    in which case it is only necessary to pass the arrays once.  In the
    more general case of arbitrary time-dependencies of mu and sigma,
    this would have to be made a bit more complicated, or the arrays
    would have to be passed separately for the two parameters.  However,
    with this idea it would not be possible anymore to take a seasonal
    oscillation with a phase.

Written by Markus Reinert, August 2020, February 2021.
"""

import calendar
from enum import Enum

import numpy as np
from scipy import stats


# Minimum value of xi (in magnitude) for which GEV is used instead of Gumbel
XI_THRESHOLD = 1e-10

# Frequency corresponding to a seasonal oscillation in the unit 1/second
SEASONAL_FREQUENCY = 2 * np.pi / (3600 * 24 * 365.2425)

# Scale factor applied to the trend for converting from cm/second to cm/century
TREND_SCALE = 3600 * 24 * 365.2425 * 100


class Modifiers(Enum):
    """Modifiers that can be used to make the GEV model time-dependent."""
    LINEAR_TREND = "linear trend"
    QUADRATIC_TREND = "quadratic trend"
    SEASONAL_OSCILLATION = "seasonal oscillation"  # or "annual cycle"
    SEASONAL_WITH_PHASE = "seasonal with phase"
    SEMIANNUAL_CYCLE = "semiannual cycle"
    CLIMATE_INDEX = "climate index"


def negative_log_likelihood(
        params, z_array, t_array=None, mu_modifiers=[], sigma_modifiers=[], fixed_xi=None,
        climate_index_array=None, verbose=False,
) -> float:
    """Calculate the negative log-likelihood of ‘z_array’ being GEV distributed.

    The parameters of the GEV distribution are given in ‘params’ in the
    order mu, sigma, xi, corresponding to location, scale, and shape
    parameter, respectively.  Alternatively, the shape can be specified
    with the parameter ‘fixed_xi’, which is especially useful in the
    context of optimisation, since the shape parameter can be difficult
    to estimate.  Using fixed_xi=0 limits the GEV family of
    distributions to the special case of a Gumbel distribution.

    The location and scale parameters of the GEV distribution can
    dependend on time and/or a climate index.  For this purpose, provide
    a ‘t_array’ corresponding to the ‘z_array’, as well as modifiers for
    mu and/or sigma (the t_array is necessary, even if there is only a
    dependence on the climate index).  In this case, the list of
    parameters must be extended by the necessary coefficients (one for a
    linear trend or a climate index, two for a seasonal oscillation) in
    the order corresponding to that of the modifiers.

    If an automatic optimisation of this function does not work, use
    verbose=True to see the params with which this function is called.

    If the GEV parameters and the z_array form an impossible
    combination, then the likelihood is zero, so the log-likelihood is
    -infinity and the negative log-likelihood is +infinity.  In this
    case, np.infty is returned.

    Since this function may return np.infty, scipy's optimize may raise
    a “RuntimeWarning: invalid value encountered in subtract”.  This
    warning can be suppressed by wrapping the call to optimize in “with
    np.errstate(invalid="ignore"):”.  This warning does not appear when
    np.infty is replaced by a very large number, e.g. 1e200.
    """
    if verbose:
        print(params)
    # Check whether it is a normal or a time-dependent GEV and get its parameters
    if t_array is None:
        if fixed_xi is None:
            mu, sigma, xi = params
        else:
            mu, sigma = params
            xi = fixed_xi
        if sigma < 0:
            return np.infty
    else:
        mu, sigma, xi = time_dependent_GEV_parameters(
            params, t_array, mu_modifiers, sigma_modifiers, fixed_xi, climate_index_array
        )
        if sigma_modifiers and any(sigma < 0):
            return np.infty
        elif not sigma_modifiers and sigma < 0:
            return np.infty
    # Calculate and return the negative log-likelihood
    if abs(xi) < XI_THRESHOLD:
        # Gumbel
        return np.sum(
            np.log(sigma) + (z_array - mu) / sigma + np.exp(-((z_array - mu) / sigma))
        )           # Synthesis of Equations 3.9 and 6.5 of Coles (2001)
    else:
        # GEV
        # Check if any point falls out of the domain of the distribution
        if all(1 + xi * (z_array - mu) / sigma > 0):  # Equation 3.8 of Coles, 2001
            return np.sum(
                np.log(sigma)
                + (1 + 1/xi) * np.log(1 + xi * (z_array - mu) / sigma)
                + (1 + xi * (z_array - mu) / sigma) ** (-1/xi)
            )                   # Equation 6.5 of Coles (2001)
        else:
            return np.infty


def time_dependent_GEV_parameters(
        params, t_array, mu_modifiers=[], sigma_modifiers=[], fixed_xi=None,
        climate_index_array=None,
) -> "3-tuple":
    """Calculate the three GEV parameters for the given time dependence.

    See negative_log_likelihood for more information on the arguments.
    """
    i = 0
    # Parse mu
    mu = params[i]
    i += 1
    for modifier in mu_modifiers:
        if modifier == Modifiers.LINEAR_TREND:
            mu += params[i] * t_array / TREND_SCALE
            i += 1
        elif modifier == Modifiers.QUADRATIC_TREND:
            mu += params[i] * (t_array / TREND_SCALE) ** 2
            i += 1
        elif modifier == Modifiers.SEASONAL_OSCILLATION:
            mu += params[i] * np.cos(t_array * SEASONAL_FREQUENCY)
            i += 1
            mu += params[i] * np.sin(t_array * SEASONAL_FREQUENCY)
            i += 1
        elif modifier == Modifiers.SEASONAL_WITH_PHASE:
            mu += params[i] * np.cos(t_array * SEASONAL_FREQUENCY + params[i+1])
            i += 2
        elif modifier == Modifiers.SEMIANNUAL_CYCLE:
            mu += params[i] * np.cos(t_array * 2 * SEASONAL_FREQUENCY)
            i += 1
            mu += params[i] * np.sin(t_array * 2 * SEASONAL_FREQUENCY)
            i += 1
        elif modifier == Modifiers.CLIMATE_INDEX:
            mu += params[i] * climate_index_array
            i += 1
        else:
            raise NotImplementedError("unknown modifier {!r} for mu".format(modifier))
    # Parse sigma
    sigma = params[i]
    i += 1
    for modifier in sigma_modifiers:
        if modifier == Modifiers.LINEAR_TREND:
            sigma += params[i] * t_array / TREND_SCALE
            i += 1
        elif modifier == Modifiers.QUADRATIC_TREND:
            sigma += params[i] * (t_array / TREND_SCALE) ** 2
            i += 1
        elif modifier == Modifiers.SEASONAL_OSCILLATION:
            sigma += params[i] * np.cos(t_array * SEASONAL_FREQUENCY)
            i += 1
            sigma += params[i] * np.sin(t_array * SEASONAL_FREQUENCY)
            i += 1
        elif modifier == Modifiers.SEASONAL_WITH_PHASE:
            sigma += params[i] * np.cos(t_array * SEASONAL_FREQUENCY + params[i+1])
            i += 2
        elif modifier == Modifiers.SEMIANNUAL_CYCLE:
            sigma += params[i] * np.cos(t_array * 2 * SEASONAL_FREQUENCY)
            i += 1
            sigma += params[i] * np.sin(t_array * 2 * SEASONAL_FREQUENCY)
            i += 1
        elif modifier == Modifiers.CLIMATE_INDEX:
            sigma += params[i] * climate_index_array
            i += 1
        else:
            raise NotImplementedError("unknown modifier {!r} for sigma".format(modifier))
    # Parse xi
    if fixed_xi is None:
        xi = params[i]
        i += 1
    else:
        xi = fixed_xi
    # Check that all parameters have been parsed
    if len(params) > i:
        raise ValueError("expected {} parameters but {} were given".format(i, len(params)))
    return mu, sigma, xi


def get_year_selection(year, time_array):
    """Get the Boolean array that selects ‘year’ from ‘time_array’.

    The time in ‘time_array’ must be in seconds since the Epoch.

    This function is useful for the selection of annual maxima from a
    time series.
    """
    t_start = calendar.timegm((year, 1, 1, 0, 0, 0, -1, -1, -1))
    t_end = calendar.timegm((year+1, 1, 1, 0, 0, 0, -1, -1, -1))
    return (time_array >= t_start) & (time_array < t_end)


def get_month_selection(year, month, time_array):
    """Get the Boolean array that selects ‘year’-‘month’ from ‘time_array’.

    The month must be between 1 and 12.
    The time in ‘time_array’ must be in seconds since the Epoch.

    This function is useful for the selection of monthly maxima from a
    time series.
    """
    t_start = calendar.timegm((year, month, 1, 0, 0, 0, -1, -1, -1))
    if month == 12:
        t_end = calendar.timegm((year+1, 1, 1, 0, 0, 0, -1, -1, -1))
    else:
        t_end = calendar.timegm((year, month+1, 1, 0, 0, 0, -1, -1, -1))
    return (time_array >= t_start) & (time_array < t_end)


def check_significance(D, k, alpha=0.05):
    """Perform a test for significance based on the deviance statistic.

    The test based on the deviance statistic as described by Coles
    (2001) compares the log-likelihoods of two statistical models, of
    which one generalises the other.  That means, the more complex model
    must have the same parameters as the simpler model, and one or
    several more.  Thus the more complex model has always a higher (or
    equal) likelihood than the simpler model.  The increase of
    log-likelihood is compared to a chi-squared distribution to check if
    the complex model significantly improves the simpler model.

    The result of the statistical test is printed as a text and returned
    as either True (increase in likelihood is significant) or False
    (more complex model is not significantly better).

    Input parameters (all non-negative):
    D:  2-times the difference of log-likelihood,
    k:  difference in the number of parameters,
    alpha:  significance level.

    """
    print(
        "Does the more complex model significantly improve the simpler model (alpha = {})?"
        .format(alpha)
    )
    D_min = stats.chi2.ppf(1 - alpha, df=k)
    if D > D_min:
        print("Yes!  D = {:.2f} > {:.2f}".format(D, D_min))
        return True
    else:
        print("No.  D = {:.2f}, but must be more than {:.2f} for significance".format(D, D_min))
        return False


def compute_amplitude_and_phase(coeff_cos, coeff_sin, std_cos, std_sin):
    """Compute amplitude and phase of a*cos(t) + b*sin(t).

    This function can be used to convert an expression of the form
    a*cos(t) + b*sin(t) to amp * cos(t - phase).

    The input parameters are coeff_cos (=a) and coeff_sin (=b) and their
    respective standard errors.

    The return value is (amp, phase, std_amp, std_phase), where:
     - amp is the amplitude sqrt(a**2 + b**2),
     - phase is such that tan(phase) = b/a,
     - std_amp is the standard error of the amplitude,
     - std_phase is the standard error of the phase.

    """
    amp = np.sqrt(coeff_cos ** 2 + coeff_sin ** 2)
    phase = -np.arctan2(-coeff_sin, coeff_cos)
    std_amp = np.sqrt(
        (std_cos * coeff_cos / amp) ** 2 + (std_sin * coeff_sin / amp) ** 2
    )
    std_phase = np.sqrt(
        (std_cos * coeff_sin / amp**2)**2 + (std_sin * coeff_cos / amp**2)**2
    )
    return amp, phase, std_amp, std_phase
