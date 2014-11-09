"""An electrophoresis solver."""
# pylint: disable=W0212
import numpy as np
import scipy.integrate as integrate
from scipy.signal import gaussian
from collections import OrderedDict


class Migrate(object):

    """A class for performing electromigration calculations."""

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
    solution_dict = OrderedDict()
    full_solution = []
    full_solution_dict = OrderedDict()
    solver_info = None
    pH = 7
    epsilon = 0.75
    NI = 10
    Kag = 0.01
    pointwave = 1e-4
    t = 0
    adaptive_grid = True
    calls = 0
    u = 0
    N_window = 10

    def __init__(self, system):
        """Initialize with a system from the constructor class."""
        self.x = np.array(system.domain)
        self.z = np.array(self.x[:])
        self.N = self.x.size
        self.set_dz()
        self.ions = system.ions
        self.M = len(self.ions)
        self.set_ion_properties()
        self.initial_concentrations = np.array(system.concentrations)
        self.concentrations = self.initial_concentrations
        self.V = system.V
        self.t = 0.0
        self.differ = self.Differentiate(self.N, self.dz, method='6th-Order')
        self.differ_dissipative = self.Differentiate(self.N, self.dz, method='dissipative')

    def first_derivative(self, x_input, mode='compact'):
        """Calculate the first derivative with respect to z."""
        if mode == 'compact':
            return self.differ.first_derivative(x_input.T).T
        else:
            return self.differ_dissipative.first_derivative(x_input.T).T

    def second_derivative(self, x_input, mode='compact'):
        """Calculate the second derivative with respect to z."""
        if mode == 'compact':
            return self.differ.second_derivative(x_input.T).T
        else:
            return self.differ_dissipative.second_derivative(x_input.T).T

    def set_ion_properties(self):
        """Set the properties of ions in the system."""
        self.diffusivity = np.array([[ion.diffusivity(self.pH)]
                                    for ion in self.ions])
        self.mobility = np.array([[ion.effective_mobility(self.pH)]
                                  for ion in self.ions])
        self.molar_conductivity = np.array([[ion.molar_conductivity(self.pH)]
                                            for ion in self.ions])

    def set_dz(self):
        """Set spatial step size in the z domain."""
        self.dz = self.z[1]-self.z[0]

    def set_current(self, concentrations):
        """Calculate the current based on a fixed voltage drop."""
        if self.adaptive_grid is True:
            self.j = self.V/sum(self.dz * self.first_derivative(self.x) /
                                self.conductivity(concentrations))
        else:
            self.j = self.V/sum(self.dz / self.conductivity(concentrations))

    def conductivity(self, concentrations):
        """Calculate the conductivty at each location."""
        conductivity = np.sum(np.tile(self.molar_conductivity,
                                      (1, self.N))
                              * concentrations, 0)
        return conductivity

    def set_E(self, concentrations):
        """Calculate the electric field at each node."""
        self.set_current(self.concentrations)
        self.E = -self.j/self.conductivity(self.concentrations)
        if self.adaptive_grid is True:
            self.E = self.E * self.first_derivative(self.x)

    def flux(self):
        """Calculate the flux of chemical species."""
        self.set_E(self.concentrations)
        if self.adaptive_grid is True:
            cD = self.tile_n(self.diffusivity) * self.concentrations
            diffusion = \
                (self.second_derivative(cD) -
                 self.first_derivative(cD) * self.tile_m(self.xzz / self.xz)) / \
                 self.tile_m(self.xz**2)

            advection = (self.tile_n(self.mobility) * self.concentrations) *\
                self.tile_m(self.first_derivative(self.E) -
                    (self.xzz/self.xz) * self.E) + self.first_derivative(
                    (self.tile_n(self.mobility) * self.concentrations)) *\
                    self.tile_m(self.E)
            advection /= self.tile_m(self.xz**2)

            node_movement = self.tile_m((self.node_flux()-self.u) / self.xz) * \
                self.first_derivative(self.concentrations)

        else:
            diffusion = \
                self.second_derivative(np.tile(self.diffusivity,
                                               (1, len(self.z)))
                                       * self.concentrations
                                       )
            advection = \
                self.first_derivative(np.tile(self.mobility,
                                               (1, len(self.z)))
                                       * self.concentrations *
                                       self.E
                                       )
            node_movement = self.first_derivative(-self.u * self.concentrations)

        total_flux = diffusion + advection + node_movement
        # total_flux[:, -1] *= 0
        # total_flux[:, 0] *= 0

        return total_flux

    def node_flux(self):
        """Calculate the flux of nodes."""
        if self.adaptive_grid is True:
            flux = self.pointwave *\
                self.first_derivative(self.node_cost() *
                                      self.first_derivative(self.x), mode=2)
            if False:
                window = gaussian(self.N, self.N*0.1)
            else:
                window = np.bartlett(self.N_window)
            flux = np.convolve(flux, window, 'same')
            flux[0, ] = flux[-1, ] = 0
        else:
            flux = np.zeros(self.x.shape)
        return flux

    def node_cost(self):
        """Calculate the cost function of each node."""
        deriv = np.abs(self.first_derivative(self.concentrations))
        cost = deriv / np.tile(np.nanmax(deriv, 1), (len(self.z), 1)).T
        cost = np.nanmax(cost, 0) + self.Kag
        return cost

    def decompose_state(self, state):
        """Decompose the state into X and concentrations."""
        x = state[:self.N]
        concentrations = state[self.N:].reshape(self.concentrations.shape)
        return (x, concentrations)

    def compose_state(self, x, concentrations):
        """Compose X and concentrations into a state."""
        x = x.flatten()
        concentrations = concentrations.flatten()
        state = np.concatenate((x, concentrations))
        return state

    def reshaped_flux(self, t, state):
        """1-D flux function for ode solver."""
        if not t == self.t:
            self.calc_equilibrium()
            self.t = t
        (self.x, self.concentrations) = self.decompose_state(state)
        self.set_derivatives()

        x_flux = self.node_flux()
        ion_flux = self.flux()
        flux = self.compose_state(x_flux, ion_flux)
        return flux

    def set_derivatives(self):
        self.xz = self.first_derivative(self.x)
        self.xzz = self.second_derivative(self.x)
        # self.phiz =
        # self.phizz =

    def tile_m(self, tileable):
        return np.tile(tileable, (self.M, 1))

    def tile_n(self, tileable):
        return np.tile(tileable, (1, self.N))

    def write_solution(self, t, state, full=True, use_dict=True):
        """Write the current state to solutions."""
        (x, concentrations) = self.decompose_state(state)
        if full is False:
            self.solution.append([t, x, concentrations])
        else:
            self.full_solution.append([t, x, concentrations])

        if use_dict:
            if full:
                if t not in self.full_solution_dict.keys():
                    self.full_solution_dict[t] = (x, concentrations)
            else:
                if t not in self.solution_dict.keys():
                    self.solution_dict[t] = (x, concentrations)
        return None

    def solve(self, tmax, dt=1, method='rk45'):
        """Solve for a series of time points using an ODE solver."""
        self.solution = []
        self.full_solution = []
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

        if solver._integrator.supports_solout:
            solver.set_solout(self.write_solution)

        self.x = self.z[:]
        solver.set_initial_value(self.compose_state(self.x, self.initial_concentrations))

        while solver.successful() and solver.t < tmax:
            tnew = solver.t + dt
            if tnew > tmax:
                tnew = tmax
            solver.integrate(tnew)
            self.write_solution(solver.t, solver.y, False)

        if not solver.successful():
            print 'solver failed at time', solver.t
        else:
            (self.x, self.concentrations) = self.decompose_state(solver.y)


    def calc_equilibrium(self):
        pass
