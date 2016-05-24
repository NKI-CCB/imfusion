# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import numpy as np

try:
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    pandas2ri.activate()

    # Note robject_translations is no longer needed in Rpy2 2.6+,
    # but is kept for compatibility with older Rpy2 versions.
    r_mass = importr('MASS')
    r_stats = importr('stats',
                      robject_translations={'format_perc': '_format_perc'})
except ImportError:
    r_mass = None
    r_stats = None


class NegativeBinomial(object):

    """ Models a negative binomial distribution.

    NegativeBinomial class that wraps functionality from the R MASS
    and stats packages for fitting and evaluating Negative Binomials.
    """

    def __init__(self, mu=None, size=None, data=None):
        if r_mass is None or r_stats is None:
            raise ValueError('Rpy2 is required to use the '
                             'NegativeBinomial test')

        if data is not None:
            mu, size = self._fit(data)

        self._mu = mu
        self._size = size

    @property
    def is_discrete(self):
        return True

    def fit(self, data):
        mu, size = self._fit(data)
        self._mu = mu
        self._size = size

    def pdf(self, x, log=False):
        self._check_if_fit()
        return np.array(r_stats.dnbinom(x=x, mu=self._mu,
                                        size=self._size, log=log))

    def cdf(self, q, log_p=False):
        self._check_if_fit()
        return self._pnbinom(q, mu=self._mu, size=self._size,
                             lower_tail=True, log_p=log_p)

    def sf(self, q, log_p=False):
        self._check_if_fit()
        return self._pnbinom(q, mu=self._mu, size=self._size,
                             lower_tail=False, log_p=log_p)

    def _check_if_fit(self):
        if self._mu is None or self._size is None:
            raise ValueError('Distribution has not yet been fit')

    @staticmethod
    def _fit(data):
        fit = r_mass.fitdistr(data, densfun='negative binomial')
        fit = dict(zip(fit.names, [np.array(item) for item in fit]))
        mu, size = fit['estimate'][1], fit['estimate'][0]

        return mu, size

    @staticmethod
    def _pnbinom(q, mu, size, **kwargs):
        p_val = r_stats.pnbinom(q, mu=mu, size=size, **kwargs)
        return p_val[0] if len(p_val) == 1 else np.array(p_val)

    def __repr__(self):
        return '<NegativeBinomial mu={}, size={}>'.format(self._size, self._mu)

    def __str__(self):
        return self.__repr__()
