"""A module for calculating the roots of a set of polynomials."""
import numpy as np
from scipy import optimize
import warnings


class Multiroot(object):
    """A class for fast roots-finding for an array of polynomials."""

    method = 'hybr'
    use_jac = True
    max_nodes = 50
    enforce_positive = True
    solver_options = {'band': (0, 0),
                      'col_deriv': True
                      }

    def __init__(self):
        pass

    def __call__(self, polys, guess=None):
        if guess is None:
            return self._analytical_solve(polys)
        elif polys.shape[1] > self.max_nodes:
            return self._split_solve(polys, guess)
        else:
            return self._solve(polys, guess)

    def _split_solve(self, polys, guess):
        # Calculate split indices
        splits = -(-polys.shape[1]//self.max_nodes)
        split_size = polys.shape[1]//splits
        split_indices = [i * split_size for i in range(1, splits)]

        # Split the domain into sections
        split_polys = np.hsplit(polys, split_indices)
        split_guess = np.hsplit(guess, split_indices)

        # Calculate the section results and cat.
        split_roots = [self._solve(sub_poly, sub_guess)
                       for sub_poly, sub_guess
                       in zip(split_polys, split_guess)
                       ]
        return np.concatenate(split_roots)

    def _solve(self, polys, guess):
        roots = self._optimize_solve(polys, guess)
        return self._ensure_positive(roots, polys)

    def _optimize_solve(self, polys, guess):
        offset = polys.shape[0]//2
        n = np.arange(polys.shape[0], 0., -1.) - offset

        objective = self._get_objective(polys, n)
        jac = self._get_jacobian(polys, n, offset)

        roots = optimize.root(objective, guess,
                              jac=jac, method=self.method,
                              options=self.solver_options)

        if not roots.success:
            warnings.warn(roots.message)

        return roots.x

    def _get_objective(self, polys, n):
        def objective(x):
            xn = x ** n[:, np.newaxis]
            return np.sum(polys * xn, 0)

        return objective

    def _get_jacobian(self, polys, n, offset):
        if not self.use_jac:
            return None
        else:
            p2 = polys[:, :] * np.arange(polys.shape[0]-offset,
                                         -offset,
                                         -1)[:, np.newaxis]
            m = np.arange(p2.shape[0], 0., -1)-offset-1

            def jac(x):
                xn = x ** m[:, np.newaxis]
                return np.diag(np.sum(p2 * xn, 0))

            return jac

    def _analytical_solve(self, polys):
        return np.apply_along_axis(self._1d_analytical_solve,
                                   axis=0, arr=polys)

    def _1d_analytical_solve(self, subpoly):
        root = np.roots(subpoly)
        root = [r for r in root if r.real > 0 and r.imag == 0]
        if root:
            if len(root) != 1:
                warnings.warn("Multiple roots found.")
            return float(root[0].real)
        else:
            raise RuntimeError("Failed to find pH.")

    def _ensure_positive(self, roots, polys):
        if self.enforce_positive:
            negative_roots = roots < 0
            if np.any(negative_roots):
                for idx, value in enumerate(negative_roots):
                    if value:
                        roots[idx] = self._1d_analytical_solve(polys[:, idx])
        return roots

if __name__ == '__main__':
    array = np.array([[1, 1, 3, -10],
                      [1, 8, 4, -12],
                      ])
    array = array.transpose()
    array = np.concatenate([array]*1000, axis=1)
    multiroot = Multiroot()
    guess = multiroot(array)
    print multiroot(array, guess-1)
