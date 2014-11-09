"""Defines a flux limiter class."""
import numpy as np


class Flux_Limiter(object):

    """Return a flux limiting function."""

    mode = "no_limiter"
    limiter_dict = dict()

    def __init__(self, mode):
        """Initialize."""
        pass

    def limit(self, flux):
        """Select a function and return a limiter."""
        r = self.get_r(flux)
        return self.limiter_dict[self.mode](r)

    def get_r(self, flux):
        """Calculate r, the input parameter for limiters."""
        return flux

    def minmod(self, r):
        """The minmod limiter function."""
        pass

    def superbee(self, r):
        """The superbee limiter function."""
        pass

    def no_limiter(self, r):
        """Return zeros, no limiter function."""
        return np.zeros(r.shape)

    limter_dict = {"no_limiter": no_limiter,
                   "minmod": minmod,
                   "superbee": superbee
                   }
