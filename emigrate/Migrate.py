import numpy as np
import scipy.integrate as integrate
from scipy.signal import gaussian


class Migrate(object):
    import constants
    from Differentiate import Differentiate
    V = None
    E = None
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
    pH = 7
    epsilon = 0.75
    Kag = 0.01
    pointwave = 1

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
        self.diffusivity = np.array([[ion.diffusivity(self.pH)]
                                    for ion in self.ions])
        self.mobility = np.array([[ion.effective_mobility(self.pH)]
                                  for ion in self.ions])
        self.molar_conductivity = np.array([[ion.molar_conductivity(self.pH)]
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

    def node_flux(self, concentrations):
        flux = -self.pointwave *\
            self.first_derivative(self.node_cost(concentrations) *
                                  self.first_derivative(self.x))
        flux = np.convolve(flux, gaussian(self.N, self.N*0.01), 'same')
        return flux

    def node_cost(self, concentrations):
        deriv = np.abs(self.first_derivative(concentrations))
        cost = deriv / np.tile(np.nanmax(deriv, 1), (len(self.z), 1)).T
        cost = np.nanmax(cost, 0) + self.Kag
        return cost

    def reshaped_flux(self, t, concentrations):
        if not t == self.t:
            self.calc_equilibrium()
        concentrations = concentrations.reshape(self.concentrations.shape)
        flux = np.ravel(self.flux(concentrations))
        return flux

    def solve(self, t, method='rk45'):
        self.solution = []
        solver = integrate.ode(self.reshaped_flux)

        if method == 'lsoda':
            solver.set_integrator('lsoda')

        elif method == 'rk45':
            solver.set_integrator('dopri5')

        elif method == 'rk8':
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
