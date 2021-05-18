"""Functions for the handling of parameter estimates with uncertainties.

Written by Markus Reinert, April 2020, May 2021.
"""

import itertools

import numpy as np


def get_error_bounds(f, x, params, errors, n_sigma=1):
    """Get the upper and lower bound for f(x, *(params ± n_sigma*errors)).

    The function f(x, *params) is evaluated for every possible
    combination of adding or substracting the errors (multiplied by
    n_sigma) to/from their corresponding parameters.  Two graphs (as
    functions of x) are returned that describe the upper and lower bound
    for the calculated graphs.

    If the result of this function is named "bounds", then the error
    area can be drawn in matplotlib with: plt.fill_between(x, *bounds)

    The lists params and errors must have the same number of elements.
    """
    assert len(params) == len(errors)
    upper_graph = lower_graph = f(x, *params)
    # Get every possible combination of positive and negative errors
    for error_combi in itertools.product(*[[e, -e] for e in errors]):
        param_variant = [p + n_sigma * e for p, e in zip(params, error_combi)]
        graph = f(x, *param_variant)
        upper_graph = np.fmax(upper_graph, graph)
        lower_graph = np.fmin(lower_graph, graph)
    return upper_graph, lower_graph


def format_GEV_parameters(parameters, errors, join_str=", "):
    """Create a string that contains the GEV parameters in a nice format."""
    names = ["\\mu", "\\sigma", "\\xi"]
    # Write the error with one significant digit
    digits = [int(-np.floor(np.log10(e))) if e < 1 else 0 for e in errors]
    return join_str.join(
        "${name} = {param:.{digits}f} \\pm {std:.{digits}f}$".format(
            name=n, param=p, std=e, digits=d
        ) for n, p, e, d in zip(names, parameters, errors, digits)
    )
