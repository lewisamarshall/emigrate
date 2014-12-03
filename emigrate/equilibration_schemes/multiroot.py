"""A function for calculating the roots of a set of polynomials."""
from numpy.polynomial.polynomial import polycompanion
from scipy.linalg import block_diag
from scipy.sparse.linalg import eigs
import numpy


def multiroot(polys, last_roots=None):
    """A function for calculating the roots of a set of polynomials."""
    comp = [polycompanion(p) for p in polys]
    n = len(comp)
    all_comp = block_diag(*comp)
    if last_roots is not None:
        w, v = eigs(all_comp, n, sigma=0, v0=last_roots, OPpart='r')
    else:
        w, v = eigs(all_comp, n, sigma=0, OPpart='r')

    return w, v

if __name__ == '__main__':
    from numpy.polynomial import Polynomial as P

    p = [[-1, 0, 1], [-1, 0, 1]]
    print p
    print multiroot(p)
    print [numpy.roots(pl) for pl in p]
