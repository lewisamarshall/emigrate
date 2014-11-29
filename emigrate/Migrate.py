"""An electrophoresis solver."""
import numpy as np
import scipy.integrate as integrate
from collections import OrderedDict
# pylint: disable=W0212


class Migrate(object):

    """A class for performing electromigration calculations."""

    import constants
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
    atol = 1e-12
    rtol = 1e-6
    equlibrator = None
    flux_calculator = None

    def __init__(self, system, flux_mode='compact', equilibrium_mode='pH'):
        """Initialize with a system from the constructor class."""
        self.x = np.array(system.domain)
        self.z = self.x.copy()
        self.N = self.x.size
        self.set_dz()
        self.ions = system.ions
        self.M = len(self.ions)
        self.initial_concentrations = np.array(system.concentrations)
        self.concentrations = self.initial_concentrations
        self.V = system.V
        self.t = 0.0
        self.equilibrum_mode = equilibrium_mode
        self.set_equilibrium_mode()
        self.flux_mode = flux_mode
        self.set_flux_mode()

    def set_equilibrium_mode(self):
        """Import an equilibration object to calculate ion properties."""
        if self.equilibrum_mode == 'fixed':
            from equilibration_schemes import Fixed
            self.equlibrator = Fixed(self.ions, self.pH, self.concentrations)
        elif self.equilibrum_mode == 'pH':
            from equilibration_schemes import Variable_pH
            self.equlibrator = Variable_pH(self.ions, self.pH,
                                           self.concentrations)
        else:
            pass
        self.mobility, self.diffusivity, self.molar_conductivity = \
            self.equlibrator.equilibrate(self.concentrations)

    def set_flux_mode(self):
        """Import a flux calculator to calculate ion fluxes."""
        if self.flux_mode == 'compact':
            from flux_schemes import Compact
            self.flux_calculator = Compact(self.N,
                                           self.dz,
                                           self.V,
                                           self.mobility,
                                           self.diffusivity,
                                           self.molar_conductivity)
        else:
            pass

    def set_dz(self):
        """Set spatial step size in the z domain."""
        self.dz = self.z[1]-self.z[0]

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

        solver = integrate.ode(self.objective)

        solver.set_integrator(method,
                              atol=self.atol,
                              rtol=self.rtol
                              )
        if solver._integrator.supports_solout:
            solver.set_solout(self.solout)

        solver.set_initial_value(self.compose_state(self.x,
                                                    self.initial_concentrations
                                                    ))

        while solver.successful() and solver.t < tmax:
            tnew = solver.t + dt
            if tnew > tmax:
                tnew = tmax
            solver.integrate(tnew)
            self.write_solution(solver.t, solver.y, False)
        (self.x, self.concentrations) = self.decompose_state(solver.y)

        if not solver.successful():
            print 'solver failed at time', solver.t

    def solout(self, t, state):
        """Perform actions when a successful solution step is found."""
        (self.x, self.concentrations) = self.decompose_state(state)
        self.write_solution(t, state)
        self.equlibrator.equilibrate(self.concentrations)

    def objective(self, t, state):
        """The objective function of the solver."""
        self.t = t
        (self.x, self.concentrations) = self.decompose_state(state)
        ion_flux = self.flux_calculator.flux(self.x, self.concentrations)
        x_flux = np.zeros(self.x.shape)
        flux = self.compose_state(x_flux, ion_flux)
        return flux
