"""Create the differentiate class to take derivatives quickly."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg


class Differentiate(object):

    """Take derivatives quickly."""

    A1 = None
    B1 = None
    A2 = None
    B2 = None
    sparse = True
    factorized = True

    def __init__(self, N, dz, method):
        """Initialize with a length and step size."""
        self.N = N
        self.dz = dz
        self.method = method
        self.set_matrices()

    def first_derivative(self, x):
        """Take the first derivative of the input."""
        if self.sparse is True:
            if self.factorized is True:
                derivative = self.fA1(self.B1.dot(x))
            else:
                derivative = linalg.spsolve(self.A1, self.B1.dot(x))
        else:
            derivative = np.linalg.solve(self.A1, np.dot(self.B1, x))
        return derivative

    def second_derivative(self, x):
        """Take the second derivative of the input."""
        if self.sparse is True:
            if self.factorized is True:
                derivative = self.fA2(self.B2.dot(x))
            else:
                derivative = linalg.spsolve(self.A2, self.B2.dot(x))
        else:
            derivative = np.linalg.solve(self.A2, np.dot(self.B2, x))
        return derivative

    def set_matrices(self):
        """Set up all required matrices."""
        self.set_A1()
        self.set_A2()
        self.set_B1()
        self.set_B2()
        if self.factorized is True:
            self.fA1 = linalg.factorized(self.A1)
            self.fA2 = linalg.factorized(self.A2)

    def set_A1(self):
        """Setup for A1."""
        if self.method == '6th-Order':
            internal_function = [1./3., 1., 1./3.]
            boundary_functions = [[1., 4.], [1./6., 1., 1./2.]]
        elif self.method == 'dissipative':
            pass

        self.A1 = self.construct_matrix(boundary_functions, internal_function)

    def set_B1(self):
        """Setup for B1."""
        if self.method == '6th-Order':
            internal_function = [-1./36., -14./18., 0., 14./18., 1./36.]
            boundary_functions = [[-37./12., 2./3., 3., -2./3., 1./12.],
                                  [-10./18., -1./2., 1., 1./18.]]
        elif self.method == 'dissipative':
            pass

        self.B1 = self.construct_matrix(boundary_functions,
                                        internal_function,
                                        True)
        self.B1 /= self.dz

    def set_A2(self):
        """Setup for A2."""
        if self.method == '6th-Order':
            internal_function = [2./11., 1., 2./11.]
            boundary_functions = [[1., 137./13.], [1./10., 1., -7./20.]]
        elif self.method == 'dissipative':
            pass

        self.A2 = self.construct_matrix(boundary_functions, internal_function)

    def set_B2(self):
        """Setup for B2."""
        if self.method == '6th-Order':
            internal_function = \
                [3./44., 12./11., -6./44.-24./11., 12./11., 3./44.]
            boundary_functions = [[1955./156., -4057./156., 1117./78., -55./78., -29./156., 7./156.],
                                  [99./80., -3., 93./40., -3./5., 3./80]]
                                  # note typo in paper saying -3/80
        elif self.method == 'dissipative':
            pass

        self.B2 = self.construct_matrix(boundary_functions,
                                        internal_function,
                                        False)
        self.B2 /= self.dz**2

    def construct_matrix(self, boundary_functions,
                         internal_function, invert=False):
        """Construct matrices based on inputs."""
        N = self.N
        l = len(internal_function)
        construct = [[0]*i + internal_function +
                     [0] * (N-i+1) for i in range(N)]
        construct = np.array(construct)[:, (l-1.)/2.:-(l+3.)/2.]
        construct[0, :] = boundary_functions[0] + \
            [0] * (N - len(boundary_functions[0]))
        construct[1, :] = boundary_functions[1] + \
            [0] * (N - len(boundary_functions[1]))
        boundary_functions[0].reverse()
        boundary_functions[1].reverse()
        construct[-1, :] = [0] * (N - len(boundary_functions[0])) + \
            boundary_functions[0]
        construct[-2, :] = [0] * (N - len(boundary_functions[1])) + \
            boundary_functions[1]
        if invert is True:
            construct[-2:, :] = -construct[-2:, :]

        if self.sparse is True:
            construct = sp.csc_matrix(construct)

        return construct

if __name__ == '__main__':
    from scipy.special import erf
    from matplotlib import pyplot as plot
    Nt = 30
    z = np.linspace(-1, 1, Nt)
    x = np.array([erf(z*3), erf(z*2)])
    my_diff = Differentiate(Nt, 1, method='6th-Order')
    # print my_diff.A1, '\n'
    # print my_diff.B1

    if False:
        plot.plot(z, x)
        plot.show()

    if True:
        d1 = my_diff.first_derivative(x.T)
        d2 = my_diff.second_derivative(x.T)

    if True:
        # plot.plot(z, np.ravel(x[0,:]))
        plot.plot(z, np.ravel(d2[:, 0]))
        plot.plot(z, np.ravel(d2[:, 1]))
        plot.plot(z, np.ravel(d1[:, 0]))
        plot.plot(z, np.ravel(d1[:, 1]))
        plot.show()
