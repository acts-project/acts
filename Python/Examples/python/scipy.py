# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""ROOT-free Gaussian fit backend for `ActsExamples::HistogramFitFunction`,
using `scipy.optimize.curve_fit`. Suitable for use in the PyPI distribution,
where the ROOT plugin (and `ActsPlugins::RootHistogramFit`) is not available.
"""

import numpy as np
from scipy.optimize import curve_fit


def _gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def _gaussian_jac(x, amplitude, mean, sigma):
    """Analytic Jacobian of `_gaussian` w.r.t. (amplitude, mean, sigma).

    Routes curve_fit's underlying MINPACK call through `_lmder` (exact
    derivatives) instead of `_lmdif` (forward-difference approximation).
    """
    z = (x - mean) / sigma
    g = np.exp(-0.5 * z * z)
    return np.stack(
        [g, amplitude * g * z / sigma, amplitude * g * z * z / sigma], axis=-1
    )


def makeScipyHistogramFitFunction(maxfev=200, ftol=1e-6, xtol=1e-6):
    """Build a Gaussian-fit `HistogramFitFunction` backed by
    `scipy.optimize.curve_fit`.

    The returned callable matches `ActsExamples::HistogramFitFunction`'s
    signature: `(hist, range) -> Optional[(mean, sigma, meanError,
    sigmaError)]`. Drops empty bins rather than weighting them at sigma=1,
    mirroring ROOT's "SQ0". The analytic Jacobian plus a bounded
    `maxfev`/`ftol`/`xtol` keeps `fitProfiles()` fast versus MINPACK's
    defaults.

    @param maxfev Maximum number of function evaluations
    @param ftol Relative error desired in the sum of squares
    @param xtol Relative error desired in the approximate solution
    """

    def fit(hist, rng):
        values = hist.histogram.values()
        edges = np.asarray(hist.histogram.axis(0).edges)
        centres = 0.5 * (edges[:-1] + edges[1:])

        if rng is not None:
            xMin, xMax = rng
            mask = (centres >= xMin) & (centres <= xMax)
            centres = centres[mask]
            values = values[mask]

        if np.count_nonzero(values) < 3 or values.sum() <= 0:
            return None

        mean0 = np.average(centres, weights=np.clip(values, 0, None))
        sigma0 = max(
            np.sqrt(
                np.average((centres - mean0) ** 2, weights=np.clip(values, 0, None))
            ),
            1e-6,
        )
        amplitude0 = values.max()

        keep = values > 0
        fitCentres = centres[keep]
        fitValues = values[keep]
        errors = np.sqrt(fitValues)

        try:
            with np.errstate(all="ignore"):
                popt, pcov = curve_fit(
                    _gaussian,
                    fitCentres,
                    fitValues,
                    p0=[amplitude0, mean0, sigma0],
                    sigma=errors,
                    absolute_sigma=True,
                    jac=_gaussian_jac,
                    maxfev=maxfev,
                    ftol=ftol,
                    xtol=xtol,
                )
        except RuntimeError:
            return None

        if not np.all(np.isfinite(pcov)):
            return None

        meanError = float(np.sqrt(pcov[1, 1]))
        sigmaError = float(np.sqrt(pcov[2, 2]))
        return (float(popt[1]), abs(float(popt[2])), meanError, sigmaError)

    return fit
