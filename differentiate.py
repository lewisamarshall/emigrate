import numpy as np


class differentiator(object):

    A1 = None
    B1 = None
    A2 = None
    B2 = None

    def __init__(self, N, dz, method):
        self.N = N
        self.dz = dz
        self.method = method

        self.set_matrices()

    def first_derivative(self, x):
        derivative = np.linalg.solve(self.A1, np.dot(self.B1, x))
        return derivative

    def second_derivative(self, x):
        derivative = np.linalg.solve(self.A2, np.dot(self.B2, x))
        return derivative

    def set_matrices(self):
        self.set_A1()
        self.set_A2()
        self.set_B1()
        self.set_B2()

    def set_A1(self):
        if self.method == '6th-Order':
            internal_function = [1./3., 1, 1./3.]
            boundary_functions = [[1., 4.], [1./6., 1, 1./2.]]
        elif self.method == 'dissipative':
            pass

        self.A1 = self.construct_matrix(boundary_functions, internal_function)

    def set_A2(self):
        if self.method == '6th-Order':
            internal_function = [2./11., 1, 2./11.]
            boundary_functions = [[1., 137./13.], [1./10., 1, -7./20.]]
        elif self.method == 'dissipative':
            pass

        self.A2 = self.construct_matrix(boundary_functions, internal_function)

    def set_B1(self):
        if self.method == '6th-Order':
            internal_function = [-1./36., -14./18., 0., 14./18., 1./36.]
            boundary_functions = [[-37./12, 2./3., 3., -2./3., 1./12.],
                                  [-10./18., -1./2., 1., 1./18.]]
        elif self.method == 'dissipative':
            pass

        self.B1 = self.construct_matrix(boundary_functions, internal_function)
        self.B1 /= self.dz

    def set_B2(self):
        if self.method == '6th-Order':
            internal_function = [3./44., 12./11., -6./44.-24./11., 12./11., 3./44.]
            boundary_functions = [[1955./156., -4057./156., 1117./78., -55./78., -29./156., 7./156.],
                                  [99./80., -3., 93./40., -3./5., 3./80]] # note typo in paper saying -3/80
            # boundary_functions=[[0],[0]]
        elif self.method == 'dissipative':
            pass

        self.B2 = self.construct_matrix(boundary_functions, internal_function)
        self.B2 /= self.dz**2

    def construct_matrix(self, boundary_functions, internal_function):
        N = self.N
        l = len(internal_function)
        construct = [[0]*i + internal_function + [0] * (N-i+1) for i in range(N)]
        construct = np.array(construct)[:, (l-1.)/2.:-(l+3.)/2.]
        construct[0, :] = boundary_functions[0] + [0] * (N - len(boundary_functions[0]))
        construct[1, :] = boundary_functions[1] + [0] * (N - len(boundary_functions[1]))
        boundary_functions[0].reverse()
        boundary_functions[1].reverse()
        construct[-1, :] = [0] * (N - len(boundary_functions[0])) + boundary_functions[0]
        construct[-2, :] = [0] * (N - len(boundary_functions[1])) + boundary_functions[1]
        return construct

if __name__ == '__main__':
    from scipy.special import erf
    from matplotlib import pyplot as plot
    N = 50
    z = np.linspace(-1, 1, N)
    x = np.array([erf(z*3), erf(z*2)])
    my_diff = differentiator(N, 1, method='6th-Order')

    if False:
        plot.plot(z, x)
        plot.show()

    if True:
        d1 = my_diff.first_derivative(x.T)
        d2 = my_diff.second_derivative(x.T)

    if True:
        plot.plot(z, np.ravel(x[0,:]))
        plot.plot(z, np.ravel(d2[:,0]))
        plot.show()
