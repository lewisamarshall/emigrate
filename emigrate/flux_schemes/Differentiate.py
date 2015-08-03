"""Create the differentiate class to take derivatives quickly."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg


class Differentiate(object):

    """Take derivatives quickly."""

    A1 = None
    fA1 = None
    B1 = None
    A2 = None
    fA2 = None
    B2 = None
    M = None
    fM = None
    epsilon = 1.

    def __init__(self, N, dz, method):
        """Initialize with a length and step size."""
        self.N = N
        self.dz = dz
        self.method = method
        self.set_matrices()

    def first_derivative(self, x):
        """Take the first derivative of the input."""
        return self.fA1(self.B1.dot(x))

    def second_derivative(self, x):
        """Take the second derivative of the input."""
        return self.fA2(self.B2.dot(x))

    def smooth(self, x):
        """Smooth x using an implicit smoothing formula."""
        return self.fM(self.A2.dot(x))

    def set_matrices(self):
        """Set up all required matrices."""
        self.set_A1()
        self.set_A2()
        self.set_B1()
        self.set_B2()
        self.set_M()
        self.fA1 = linalg.factorized(self.A1)
        self.fA2 = linalg.factorized(self.A2)
        self.fM = linalg.factorized(self.M)
        self.fA2t = linalg.factorized(self.A2[1:-1, 1:-1])

    def set_A1(self):
        """Setup for A1."""
        if self.method == '6th-Order':
            internal_function = [1./3., 1., 1./3.]
            boundary_functions = [[1., 4.], [1./6., 1., 1./2.]]
            self.A1 = self.construct_matrix(boundary_functions,
                                            internal_function)
        elif self.method == 'dissipative':
            self.A1 = sp.csc_matrix(np.identity(self.N))

    def set_B1(self):
        """Setup for B1."""
        if self.method == '6th-Order':
            internal_function = [-1./36., -14./18., 0., 14./18., 1./36.]
            boundary_functions = [[-37./12., 2./3., 3., -2./3., 1./12.],
                                  [-10./18., -1./2., 1., 1./18.]]
            flag = True
        elif self.method == 'dissipative':
            internal_function = [-1./2., 0, 1./2.]
            boundary_functions = [[-1, 1]]
            flag = True

        self.B1 = self.construct_matrix(boundary_functions,
                                        internal_function,
                                        flag)
        self.B1 /= self.dz

    def set_A2(self):
        """Setup for A2."""
        if self.method == '6th-Order':
            internal_function = [2./11., 1., 2./11.]
            boundary_functions = [[1., 137./13.], [1./10., 1., -7./20.]]
            self.A2 = self.construct_matrix(boundary_functions,
                                            internal_function)
        elif self.method == 'dissipative':
            self.A2 = sp.csc_matrix(np.identity(self.N))

    def set_B2(self):
        """Setup for B2."""
        if self.method == '6th-Order':
            internal_function = \
                [3./44., 12./11., -6./44.-24./11., 12./11., 3./44.]
            boundary_functions = [[1955./156., -4057./156., 1117./78.,
                                   -55./78., -29./156., 7./156.],
                                  [99./80., -3., 93./40., -3./5., 3./80]]
            # note typo in paper saying -3/80
        elif self.method == 'dissipative':
            internal_function = [1, -2, 1]
            boundary_functions = [[1, -2, 1]]

        self.B2 = self.construct_matrix(boundary_functions,
                                        internal_function,
                                        False)
        self.B2 /= self.dz**2

    def set_M(self):
        """Set up the implicit smoothing matrix."""
        self.M = sp.lil_matrix(self.A2 - self.epsilon * self.B2)
        self.M[:, 0] = self.M[:, -1] = 0.
        self.M[0, 0] = self.M[-1, -1] = 1.
        self.M = sp.csc_matrix(self.M)

    def construct_matrix(self, boundary_functions,
                         internal_function, invert=False):
        """Construct matrices based on inputs."""
        N = self.N
        l = len(internal_function)
        construct = [[0]*i + internal_function +
                     [0] * (N-i+1) for i in range(N)]
        construct = np.array(construct)[:, (l-1.)/2.:-(l+3.)/2.]

        for idx, func in enumerate(boundary_functions):
            construct[idx, :] = func + [0] * (N - len(func))
            func.reverse()
            if invert is True:
                func = [-i for i in func]
            construct[-1-idx, :] = [0] * (N - len(func)) + func

        construct = sp.csc_matrix(construct)

        return construct

if __name__ == '__main__':
    from scipy.special import erf
    from matplotlib import pyplot as plot
    Nt = 50
    z = np.linspace(-1, 1, Nt)
    test_functions = np.array([19*erf(z*3), 25*erf(z*2), 15*(erf(z*2) +
                               .3*np.random.random(z.shape))/(10*z**2+1)])
    my_diff = Differentiate(Nt, 1, method='6th-Order')
    # my_diff = Differentiate(Nt, 1, method='dissipative')
    # print my_diff.A1.todense(), '\n'
    # print my_diff.B1.todense()

    if False:
        plot.plot(z, test_functions)
        plot.show()

    if False:
        d1 = my_diff.first_derivative(test_functions.T)
        d2 = my_diff.second_derivative(test_functions.T)

    if True:
        print my_diff.fA2
        print my_diff.fM
        print test_functions[2, :].shape
        smoothed = my_diff.smooth(test_functions[2, :])
        print smoothed.shape

    if True:
        for i in range(3):
            plot.plot(z, np.ravel(test_functions[i, :]), label='test-function')
        # plot.plot(z, np.ravel(d2[:, 0]))
        plot.plot(z, np.ravel(smoothed), label='smoothed')
        # plot.plot(z, np.ravel(d2[:, 1]))
        # plot.plot(z, np.ravel(d1[:, 0]))
        # plot.plot(z, np.ravel(d1[:, 1]))
        plot.legend(loc="upper left")
        plot.show()
