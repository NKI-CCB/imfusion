# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import numpy as np


def normalize(counts):
    """Normalizes counts for sequencing depth using the median-of-ratios."""
    size_factors = estimate_size_factors(counts)
    return counts / size_factors


def estimate_size_factors(counts):
    """Calculates size factors using the median-of-ratios approach."""

    # Convert to float.
    counts = counts.astype(float)

    # Calculate size factors.
    log_geo_means = np.mean(np.log(counts), axis=1)
    sf = np.apply_along_axis(_estimate_size_factors_col, axis=0,
                             arr=counts, log_geo_means=log_geo_means)

    return sf


def _estimate_size_factors_col(counts, log_geo_means):
    log_counts = np.log(counts)
    mask = np.isfinite(log_geo_means) & (counts > 0)
    return np.exp(np.median((log_counts - log_geo_means)[mask]))
