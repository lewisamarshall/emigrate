"""A function for calculating the roots of a set of polynomials."""
import numpy as np
from scipy.optimize import root


def multiroot(polys, last_roots, method='broyden1', use_jac=False):
    """A function for calculating the roots of a set of polynomials."""
    n = np.arange(polys.shape[0])

    if use_jac:
        p2 = polys[:, 1:] * np.arange(1, polys.shape[1])
        m = range(p2.shape[1])
        def jac(x):
            xn = x[:, np.newaxis] ** m
            return np.diag(np.sum(p2 * xn, 1))
    else:
        jac = False

    def polyval(x):
        xn = x ** n[:, np.newaxis]
        return np.sum(polys * xn, 0)

    return root(polyval, last_roots, jac=jac, method=method).x
