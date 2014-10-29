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
        derivative = np.linalg.solve(self.A1, np.dot(x_input, self.B1))
        return derivative

    def second_derivative(self, x):
        derivative = np.linalg.solve(self.A2, np.dot(x_input, self.B2))
        return derivative

    def set_matrices(self):
        self.set_A1()
        self.set_A2()
        self.set_B1()
        self.set_B2()

    def set_A1(self):
        internal_function = [1./3., 1, 1./3.]
        boundary_functions = [[1., 4.], [1./6., 1, 1./2.]]
        self.A1 = self.construct_matrix(boundary_functions, internal_function)

    def set_A2(self):
        internal_function = [2./11., 1, 2./11.]
        boundary_functions = [[1., 137./13.], [1./6., 1, 1./2.]]
        self.A2 = self.construct_matrix(boundary_functions, internal_function)


    def set_B1(self):
        internal_function = [-1./36., -14./18., 0., 14./18., 1./36.]
        boundary_functions = [[-37./12, 2./3., 3., -2./3., 1./12.],
                              [-10./18., -1./2., 1., 1./18.]]
        self.B1 = self.construct_matrix(boundary_functions, internal_function)/self.dz

    def set_B2(self):
        #internal_function = [1./3., 1, 1./3.]
        #boundary_functions = [[1., 4.], [1./6., 1, 1./2.]]
        self.B2 = self.construct_matrix(boundary_functions, internal_function)/self.dz


    def construct_matrix(self, boundary_functions, internal_function):


if __name__ == '__main__':
    from scipy.special import erf
    from matplotlib import pyplot as plot
    z = np.linspace(-1, 1, 10)
    x = erf(z*3)
    my_diff = differentiator(10, 1, method='dissipative')

    if False:
        plot.plot(z, x)
        plot.show()
