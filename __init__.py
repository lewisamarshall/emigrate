import numpy as np
import ionize
import scipy.integrate as integrate
import scipy.sparse as sparse


class Electrophoresis(object):
    import constants
    V = 100
    E = 1
    dz = None
    z = None
    x = None
    A = 0
    B = 0
    diffusivity = None
    mobility = None
    molar_conductivity = None
    j = 0
    solution = None
    solver_info = None

    def __init__(self, domain, ions, concentrations):
        self.x = np.array(domain)
        self.z = np.array(self.x[:])
        self.set_dz()
        self.ions = ions
        self.set_ion_properties()
        self.set_derivative_matrices()
        self.concentrations = np.array(concentrations)
        self.t = 0

    def first_derivative(self, x_input, method='6th-Order'):
        if method is None:
            derivative = x_input

        elif method == 'dissipative':
            derivative = []

            derivative = np.pad(np.diff(x_input, n=1, axis=1),
                                ((0, 0), (0, 1)), 'reflect') / \
                np.tile(self.dz, (len(self.ions), 1))/2 -\
                np.pad(np.diff(x_input, n=2, axis=1),
                       ((0, 0), (1, 1)), 'reflect')/2 / \
                np.tile(self.dz, (len(self.ions), 1))

        elif method == '6th-Order':
            derivative = np.linalg.solve(self.A, np.dot(self.B, np.atleast_2d(x_input).T)).T
            # derivative = np.ravel(derivative)


        return derivative

    def second_derivative(self, x_input, method='dissipative'):
        if method is None:
            derivative = input
        elif method == 'dissipative':
            derivative = self.first_derivative(self.first_derivative(x_input, 'dissipative'), 'dissipative')
        elif method == '6th-Order':
            pass
        return derivative

    def set_ion_properties(self):
        pH = 7
        self.diffusivity = np.array([[ion.diffusivity(pH)]
                                    for ion in self.ions])
        self.mobility = np.array([[ion.effective_mobility(pH)]
                                  for ion in self.ions])
        self.molar_conductivity = np.array([[ion.molar_conductivity(pH)]
                                            for ion in self.ions])

    def set_dz(self):
        self.dz = np.pad(np.diff(self.z), (0, 1), 'reflect')

    def set_current(self, concentrations):
        self.j = self.V/sum(self.dz/self.conductivity(concentrations))

    def conductivity(self, concentrations):
        conductivity = np.sum(np.tile(self.molar_conductivity,
                                      (1, len(self.dz)))
                              * concentrations, 0)
        return conductivity

    def set_E(self, concentrations):
        self.set_current(concentrations)
        self.E = self.j/self.conductivity(concentrations)

    def flux(self, concentrations):
        self.set_E(concentrations)
        diffusion = \
            self.second_derivative(np.tile(self.diffusivity,
                                           (1, len(self.z)))
                                   * concentrations
                                   )
        advection = \
            -self.first_derivative(np.tile(self.mobility,
                                           (1, len(self.z)))
                                   * concentrations *
                                   self.E
                                   )

        total_flux = diffusion + advection

        print total_flux.shape

        return total_flux

    def reshaped_flux(self, concentrations, t):
        if not t == self.t:
            self.calc_equilibrium()
        concentrations = concentrations.reshape(self.concentrations.shape)
        flux = np.ravel(self.flux(concentrations))
        return flux

    def solve(self, t):
        self.solution, self.solver_info =\
            integrate.odeint(self.reshaped_flux,
                             self.concentrations.flatten(),
                             t,
                             full_output=True)
        self.solution = [sol.reshape(self.concentrations.shape)
                         for sol in self.solution]

    def calc_equilibrium(self):
        pass

    def set_derivative_matrices(self):
        h = self.dz[0]
        N = len(self.z)
        aI = 1./3.

        A_vector = ([1./6.]*1+[1./3.]*(N-3)+[4.]*1)
        B_vectors = []
        B_vectors.append([0]*(N-5) + [1./12.])  # diag -4
        B_vectors.append([0]*(N-4) + [-2./3.])  # diag -3
        B_vectors.append([(aI*4.-1.)/12.] *(N-4) + [1./18.] + [3.])  # diag -2
        B_vectors.append([-10./18.] + [3.*aI+6. ]*(N-4) + [1.]+[2./3.])  # diag -1
        B_vectors.append([-35./12.]+[-1./2.]+[0]*(N-4)+[-1./2.]+[-35./12.])  # diag 0

        A_constructor = [A_vector+[0], np.ones(self.z.size), [0]+A_vector[::-1]]
        B_constructor = 1/h*np.array([B_vectors[0]+[0]*4,
                                      B_vectors[1]+[0]*3,
                                      B_vectors[2]+[0]*2,
                                      B_vectors[3]+[0]*1,
                                      B_vectors[4]+[0]*0,
                                      [0]*1+B_vectors[3][::-1],
                                      [0]*2+B_vectors[2][::-1],
                                      [0]*3+B_vectors[1][::-1],
                                      [0]*4+B_vectors[0][::-1],
                                      ])
        self.A = sparse.spdiags(A_constructor, range(-1, 1+1), N, N).todense()
        self.B = sparse.spdiags(B_constructor, range(-4, 4+1), N, N).todense()
        self.B = self.B / np.max(self.B*-39./12.)

if __name__ == '__main__':
    from scipy.special import erf
    from matplotlib import pyplot as plot
    my_ions = ionize.Solution(['tris', 'hydrochloric acid', 'caproic acid'],
                              [0, 0, 0]
                              ).ions[0:1]

    domain_length = 0.1
    interface_length = 0.01
    nodes = 50
    my_domain = np.linspace(-domain_length/2., domain_length/2., nodes)
    my_concentrations = np.array([np.ones(my_domain.shape)*.1,
                                 0.05-0.05*erf(my_domain/interface_length),
                                 0.05*erf(my_domain/interface_length)+0.05])
    # my_concentrations = np.array(0.05-0.05*erf(my_domain/interface_length), order=3)
    my_elec = Electrophoresis(my_domain, my_ions, my_concentrations)
    # print my_elec.A
    # print my_elec.B
    # print my_elec.concentrations.shape
    # print np.linalg.solve(my_elec.A, np.dot(my_elec.B, np.atleast_2d(my_elec.concentrations).T))
    # print my_elec.first_derivative(my_elec.concentrations)
    my_elec.solve(np.array(np.linspace(0, 2e2, 10)))
    for my_sol in my_elec.solution:
        for sub_sol in my_sol:
            # sub_sol = my_sol
            plot.plot(my_elec.z, sub_sol)
    plot.show()
