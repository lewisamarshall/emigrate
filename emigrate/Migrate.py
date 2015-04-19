"""An electrophoresis solver."""
import numpy as np
import scipy.integrate as integrate
from collections import OrderedDict
from .Electrolyte import Electrolyte
import warnings
# pylint: disable=W0212


class Migrate(object):

    """A class for performing electromigration calculations."""

    import constants

    # Migrate State
    t = 0
    V = None
    E = None
    dz = None
    z = None
    x = None
    j = 0
    pH = None
    ions = None
    concntrations = None

    # Solver Parameters
    atol = 1e-12
    rtol = 1e-6

    # Modules
    equlibrator = None
    flux_calculator = None
    adaptive_grid = False

    # Solutions
    solution = OrderedDict()
    full_solution = OrderedDict()

    def __init__(self, system, flux_mode='compact', equilibrium_mode='pH'):
        """Initialize with a system from the constructor class."""
        self.x = np.array(system.nodes)
        self.z = self.x.copy()
        self.N = self.x.size
        self.set_dz()
        self.ions = system.ions
        self.M = len(self.ions)
        self.initial_concentrations = np.array(system.concentrations)
        self.concentrations = self.initial_concentrations
        self.V = system.voltage
        self.t = 0.
        self.equilibrum_mode = equilibrium_mode
        self.set_equilibrium_mode()
        self.flux_mode = flux_mode
        self.set_flux_mode()

    def set_equilibrium_mode(self):
        """Import an equilibration object to calculate ion properties."""
        if self.equilibrum_mode == 'fixed':
            from equilibration_schemes import Fixed
            self.equlibrator_class = Fixed
        elif self.equilibrum_mode == 'pH':
            from equilibration_schemes import Variable_pH
            self.equlibrator_class = Variable_pH
        else:
            raise RuntimeError('Available equlibibrators are "fixed" and "pH".')

        self.equlibrator = self.equlibrator_class(self.ions, self.pH,
                                                  self.concentrations)
        self.equlibrator.equilibrate(self.concentrations)

    def set_flux_mode(self):
        """Import a flux calculator to calculate ion fluxes."""
        if self.flux_mode == 'compact':
            from flux_schemes import Compact
            self.flux_calculator_class = Compact
        elif self.flux_mode == 'compact adaptive':
            from flux_schemes import CompactAdaptive
            self.flux_calculator_class = CompactAdaptive
        elif self.flux_mode == 'slip':
            from flux_schemes import SLIP
            self.flux_calculator_class = SLIP
        else:
            raise RuntimeError
        self.flux_calculator = self.flux_calculator_class(self.N,
                                                          self.dz,
                                                          self.V)
        self.flux_calculator.update_ion_parameters(self.equlibrator)

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
        pH = self.equlibrator.pH
        ionic_strength = self.equlibrator.ionic_strength
        current_electrolyte = Electrolyte(nodes=x, ions=self.ions,
                                          concentrations=concentrations,
                                          pH=pH, ionic_strength=ionic_strength,
                                          voltage=self.V,
                                          current_density=self.j)
        if full:
            if t not in self.full_solution.keys():
                self.full_solution[t] = current_electrolyte
        else:
            if t not in self.solution.keys():
                self.solution[t] = current_electrolyte

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
        else:
            warnings.warn("Solver doesn't support solout.")

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
        self.equlibrator.equilibrate(self.concentrations)
        self.flux_calculator.update_ion_parameters(self.equlibrator)
        self.write_solution(t, state)

    def objective(self, t, state):
        """The objective function of the solver."""
        self.t = t
        (self.x, self.concentrations) = self.decompose_state(state)
        ion_flux = self.flux_calculator.dcdt(self.x, self.concentrations)
        x_flux = self.flux_calculator.node_flux()
        flux = self.compose_state(x_flux, ion_flux)
        return flux
