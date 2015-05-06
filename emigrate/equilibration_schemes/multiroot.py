"""A function for calculating the roots of a set of polynomials."""
import numpy as np
from scipy.optimize import root
import warnings


def multiroot(polys, last_roots, method='hybr', use_jac=True):
    """A function for calculating the roots of a set of polynomials."""

    # Choose the offset to be half the length of the polynomial
    # This helps deal with floating point overflows
    offset = polys.shape[0]//2
    # offset = 0

    # Find the length of the polynomials
    n = np.arange(polys.shape[0])-offset
    # Reverse the polynomials
    polys = polys[::-1, :]

    # Calculate the jacobians
    if use_jac:
        p2 = polys[:, :] * np.arange(-offset, polys.shape[0]-offset)[:, np.newaxis]
        m = np.arange(p2.shape[0])-offset-1
        # print p2.shape, m.shape, last_roots.shape

        def jac(x):
            xn = x ** m[:, np.newaxis]
            return np.diag(np.sum(p2 * xn, 0))

    else:
        jac = None

    def polyval(x):
        xn = x ** n[:, np.newaxis]
        return np.sum(polys * xn, 0)

    new_roots = root(polyval, last_roots, jac=jac, method=method)

    if not new_roots.success:
        warnings.warn(new_roots.message)
        # raise RuntimeError('Root finder failed. ')

    return new_roots.x
