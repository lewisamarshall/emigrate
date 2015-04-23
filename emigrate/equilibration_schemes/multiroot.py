"""A function for calculating the roots of a set of polynomials."""
import numpy as np
from scipy.optimize import root


def multiroot(polys, last_roots, method='hybr', use_jac=True):
    """A function for calculating the roots of a set of polynomials."""
    # Find the length of the polynomials
    n = np.arange(polys.shape[0])
    # Reverse the polynomials
    polys = polys[::-1, :]

    # Calculate the jacobians
    if use_jac:
        p2 = polys[1:, :] * np.arange(1, polys.shape[0])[:, np.newaxis]
        m = np.arange(p2.shape[0])

        def jac(x):
            xn = x ** m[:, np.newaxis]
            return np.diag(np.sum(p2 * xn, 0))

    else:
        jac = None

    def polyval(x):
        xn = x ** n[:, np.newaxis]
        return np.sum(polys * xn, 0)

    return root(polyval, last_roots, jac=jac, method=method).x
