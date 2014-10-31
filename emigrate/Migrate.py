import numpy as np
import ionize
import scipy.integrate as integrate


class Migrate(object):
    import constants
    from Differentiate import Differentiate
    V = .1
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
    solution = []
    solver_info = None

    def __init__(self, system):
        self.x = np.array(system.domain)
        self.z = np.array(self.x[:])
        self.N = self.x.size
        self.set_dz()
        self.ions = system.ions
        self.set_ion_properties()
        self.concentrations = np.array(system.concentrations)
        self.V = system.V
        self.t = 0.0
        self.differ = self.Differentiate(self.N, self.dz, method='6th-Order')

    def first_derivative(self, x_input):
        return self.differ.first_derivative(x_input.T).T

    def second_derivative(self, x_input):
        return self.differ.second_derivative(x_input.T).T

    def set_ion_properties(self):
        pH = 7
        self.diffusivity = np.array([[ion.diffusivity(pH)]
                                    for ion in self.ions])
        self.mobility = np.array([[ion.effective_mobility(pH)]
                                  for ion in self.ions])
        self.molar_conductivity = np.array([[ion.molar_conductivity(pH)]
                                            for ion in self.ions])

    def set_dz(self):
        self.dz = self.z[1]-self.z[0]

    def set_current(self, concentrations):
        self.j = self.V/sum(self.dz/self.conductivity(concentrations))

    def conductivity(self, concentrations):
        conductivity = np.sum(np.tile(self.molar_conductivity,
                                      (1, self.N))
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
        return total_flux

    def node_flux(self):
        pass

    def reshaped_flux(self, concentrations, t):
        if not t == self.t:
            self.calc_equilibrium()
        concentrations = concentrations.reshape(self.concentrations.shape)
        flux = np.ravel(self.flux(concentrations))
        return flux

    def ode_reshaped_flux(self, t, concentrations):
        return self.reshaped_flux(concentrations, t)

    def solve(self, t, method='rk45'):
        self.solution = []
        solver = integrate.ode(self.ode_reshaped_flux)

        if method == 'lsoda':
            solver.set_integrator('lsoda')

        elif method == 'rk45':
            solver.set_integrator('dopri5')

        elif method == 'rk8(53)':
            solver.set_integrator('dop853')

        elif method == 'vode':
            solver.set_integrator('vode')

        elif method == 'zvode':
            solver.set_integrator('zvode')

        solver.set_initial_value(self.concentrations.flatten())
        for tp in t[1:-1]:
            solver.integrate(tp)
            self.solution.append(solver.y)
            if not solver.successful():
                print 'solver failed'
                break

        self.solution = [sol.reshape(self.concentrations.shape)
                         for sol in self.solution]

    def calc_equilibrium(self):
        pass

if __name__ == '__main__':
    from scipy.special import erf
    from matplotlib import pyplot as plot
    my_ions = ionize.Solution(['tris', 'caproic acid', 'hydrochloric acid'],
                              [0, 0, 0]
                              ).ions

    domain_length = 0.1
    interface_length = 0.01
    nodes = 200
    my_domain = np.linspace(-domain_length/2., domain_length/2., nodes)
    my_concentrations = np.array([np.ones(my_domain.shape)*.1,
                                 0.05-0.05*erf(my_domain/interface_length),
                                 0.05*erf(my_domain/interface_length)+0.05])
    # my_concentrations = np.array(0.05-0.05*erf(my_domain/interface_length), order=3)
    my_elec = Migrate(my_domain, my_ions, my_concentrations)
    # print my_elec.A
    # print '\n'
    # print my_elec.B[0]
    # print my_elec.concentrations.shape
    # print np.linalg.solve(my_elec.A, np.dot(my_elec.B, np.atleast_2d(my_elec.concentrations).T))
    # print my_elec.first_derivative(my_elec.concentrations, method ='dissipative')
    deriv =  my_elec.first_derivative(my_elec.concentrations)[:,1]
    deriv = np.ravel(deriv)
    print deriv.shape, my_elec.z.shape
    # plot.plot(my_elec.z, deriv)
    # plot.plot(my_elec.z, my_elec.concentrations[1,:])
    my_elec.solve(np.array(np.linspace(0, 5e2, 10)))
    for my_sol in my_elec.solution:
        for sub_sol in my_sol:
            # sub_sol = my_sol
            # # my_elec.set_E(my_sol)
            # # print sub_sol.shape
            # deriv =  np.ravel(my_elec.first_derivative(sub_sol)[2,:])
            # conc = np.ravel(sub_sol[1,:])
            # plot.plot(my_elec.z, deriv)
            plot.plot(my_elec.z, sub_sol)
    plot.show()
#
