"""A function for calculating the roots of a set of polynomials."""
import numpy as np
from scipy.optimize import root
import warnings


def multiroot(polys, last_roots, method='hybr', use_jac=True,
              enforce_positive=True):
    """A function for calculating the roots of a set of polynomials."""

    # Choose the offset to be half the length of the polynomial
    # This helps deal with floating point overflows
    offset = polys.shape[0]//2

    # Find the length of the polynomials
    n = np.arange(polys.shape[0])-offset
    # Reverse the polynomials
    polys = polys[::-1, :]

    # Calculate the jacobians
    if use_jac:
        p2 = polys[:, :] * np.arange(-offset, polys.shape[0]-offset)[:, np.newaxis]
        m = np.arange(p2.shape[0])-offset-1

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

    x = new_roots.x

    if enforce_positive is True:
        negative_roots = x<0
        if np.any(negative_roots):
            for idx, value in enumerate(negative_roots):
                if value:
                    x[idx] = real_positive_root(polys[::-1, idx])

    return x


def real_positive_root(poly):
    cH = np.roots(poly)

    cH = [c for c in cH if c.real > 0 and c.imag == 0]
    if cH:
        if len(cH) != 1:
            warnings.warn("Multiple roots found.")
        return float(cH[0].real)
    else:
        raise RuntimeError("Failed to find pH.")
