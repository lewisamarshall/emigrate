"""An electrophoresis solver."""
# pylint: disable=W0212
import numpy as np
import scipy.integrate as integrate
# from scipy.signal import gaussian
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
    diffusivity = None
    mobility = None
    molar_conductivity = None
    j = 0
    solution = OrderedDict()
    full_solution = OrderedDict()
    solver_info = None
    pH = 7
    epsilon = 0.75
    NI = 10
    Kag = 0.01
    pointwave = 1e-5
    t = 0
    adaptive_grid = True
    calls = 0
    u = 0
    N_window = 20
    Vthermal = .025
    alpha = None
    characteristic = None
    boundary_mode = 'fixed'
    xz = None
    xzz = None

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
        self.set_Kag()

    def first_derivative(self, x_input):
        """Calculate the first derivative with respect to z."""
        return self.differ.first_derivative(x_input.T).T

    def second_derivative(self, x_input):
        """Calculate the second derivative with respect to z."""
        return self.differ.second_derivative(x_input.T).T

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

    def set_Kag(self):
        """Set the Kag parameter for spacing of low-gradient grid points."""
        self.Kag = ((self.N-self.NI)/self.NI) * self.Vthermal / self.V

    def set_alpha(self):
        """Set alpha for dissipation."""
        self.set_characteristic()
        self.alpha = 0.5 * np.maximum(np.fabs(self.characteristic/self.tile_m(self.xz)),0)

    def set_characteristic(self):
        """Calculate the characteristic speed of paramters."""
        self.characteristic = self.u + self.tile_m(self.E)*self.tile_n(self.mobility) -\
            self.node_flux()

    def set_current(self):
        """Calculate the current based on a fixed voltage drop."""
        if self.adaptive_grid is True:
            self.j = self.V/sum(self.dz * self.first_derivative(self.x) /
                                self.conductivity())
        else:
            self.j = self.V/sum(self.dz / self.conductivity())

    def conductivity(self):
        """Calculate the conductivty at each location."""
        conductivity = np.sum(np.tile(self.molar_conductivity,
                                      (1, self.N))
                              * self.concentrations, 0)
        return conductivity

    def set_E(self):
        """Calculate the electric field at each node."""
        self.set_current()
        self.E = -self.j/self.conductivity()
        if self.adaptive_grid is True:
            self.E = self.E * self.first_derivative(self.x)

    def flux(self):
        """Calculate the flux of chemical species."""
        self.set_E()
        total_flux = self.diffusive_flux() + \
            self.electromigration_flux() + \
            self.node_movement_flux()
        total_flux = self.set_boundary(total_flux)
        return total_flux

    def diffusive_flux(self):
        """Calculate flux due to diffusion."""
        cD = self.tile_n(self.diffusivity) * self.concentrations
        if self.adaptive_grid is True:
            diffusion = \
                (self.second_derivative(cD) -
                 self.first_derivative(cD) * self.tile_m(self.xzz / self.xz)) / \
                self.tile_m(self.xz**2.)
        else:
            diffusion = \
                self.second_derivative(cD)
        return diffusion

    def electromigration_flux(self):
        """Calculate flux due to electromigration."""
        uc = self.tile_n(self.mobility) * self.concentrations

        temp  = (uc * (self.first_derivative(self.E) -(self.xzz/self.xz) * self.E)
                + self.first_derivative(uc) * self.E)/self.xz**2.
        if self.adaptive_grid is True:
            electromigration = temp
        else:
            electromigration = \
                self.first_derivative(uc * self.E)
            # print np.max(electromigration-temp)/np.max(electromigration)
        return electromigration

    def node_movement_flux(self):
        if self.adaptive_grid is True:
            node_movement = self.tile_m((self.node_flux()-self.u) / self.xz) * \
                self.first_derivative(self.concentrations)
        else:
            node_movement = self.first_derivative(-self.u * self.concentrations)
        return node_movement

    def set_boundary(self, flux):
        if self.boundary_mode == 'fixed':
            flux[:, 0] *= 0
            flux[:, -1] *= 0
        elif self.boundary_mode == 'characteristic':
            pass
        return flux

    def node_flux(self):
        """Calculate the flux of nodes."""
        if self.adaptive_grid is True:
            flux = self.pointwave *\
                self.first_derivative(self.node_cost() *
                                      self.first_derivative(self.x))
            window = np.bartlett(self.N_window)
            flux = np.convolve(flux, window, 'same')
            flux[0, ] = flux[-1, ] = 0.
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

    def tile_m(self, tileable):
        return np.tile(tileable, (self.M, 1))

    def tile_n(self, tileable):
        return np.tile(tileable, (1, self.N))

    def write_solution(self, t, state, full=True):
        """Write the current state to solutions."""
        (x, concentrations) = self.decompose_state(state)

        if full:
            if t not in self.full_solution.keys():
                self.full_solution[t] = (x, concentrations)
        else:
            if t not in self.solution.keys():
                self.solution[t] = (x, concentrations)

    def solve(self, tmax, dt=1, method='dopri5'):
        """Solve for a series of time points using an ODE solver."""
        self.solution = OrderedDict()
        self.full_solution = OrderedDict()
        self.x = self.z[:]

        solver = integrate.ode(self.reshaped_flux)

        solver.set_integrator(method)
        if solver._integrator.supports_solout:
            solver.set_solout(self.solout)

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

    def solout(self, t, state):
        self.write_solution(t, state)
        self.calc_equilibrium()

    def calc_equilibrium(self):
        pass
