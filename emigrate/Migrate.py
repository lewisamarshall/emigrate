"""An electrophoresis solver."""
import numpy as np
import scipy.integrate as integrate
from collections import OrderedDict
from .Electrolyte import Electrolyte
import warnings
import flux_schemes
import equilibration_schemes
from Electromigration import Electromigration
# pylint: disable=W0212


class Migrate(object):

    """A class for performing electromigration calculations."""

    import constants

    # Migrate State
    t = 0.
    V = None
    E = None
    dz = None
    z = None
    x = None
    j = 0
    pH = None
    ions = None
    concentrations = None
    u = 0.

    # Solver Parameters
    atol = 1e-12
    rtol = 1e-6

    # Modules
    equlibrator = None
    flux_calculator = None
    adaptive_grid = False

    # Solutions
    electromigration = None

    #Frame of reference info
    frame = None
    edge = 'right'

    def __init__(self, system, flux_mode='compact', equilibrium_mode='pH'):
        """Initialize with a system from the constructor class."""

        # Prepare System Domain
        self._prep_domain(system.nodes)

        # Prepare ions
        self.ions = system.ions
        self.concentrations = np.array(system.concentrations)

        # Set system voltage mode
        self.V = system.voltage
        self.u = system.u

        # Set equilibrium mode.
        self.equilibrum_mode = equilibrium_mode
        self._set_equilibrium_mode()

        # Set flux mode
        self.flux_mode = flux_mode
        self._set_flux_mode()

        # Create empty solution dictionaries
        self.electromigration = Electromigration(self.ions)
        self._write_solution(0,
                             self._compose_state(self.x,
                                                 self.concentrations),
                             full=False)

    def _prep_domain(self, nodes):
        self.x = np.array(nodes)
        self.z = np.linspace(min(self.x), max(self.x), len(self.x))
        self.N = self.x.size
        self.dz = self.z[1]-self.z[0]

    def _set_equilibrium_mode(self):
        """Import an equilibration object to calculate ion properties."""
        if self.equilibrum_mode == 'fixed':
            self.equlibrator = equilibration_schemes.Fixed
        elif self.equilibrum_mode == 'pH':
            self.equlibrator = equilibration_schemes.Variable_pH
        else:
            raise RuntimeError('Available equlibibrators are "fixed" and "pH".'
                               )

        self.equlibrator = self.equlibrator(self.ions, self.pH,
                                            self.concentrations)
        self.equlibrator.equilibrate(self.concentrations)

    def _set_flux_mode(self):
        """Import a flux calculator to calculate ion fluxes."""
        if self.flux_mode == 'compact':
            self.flux_calculator = flux_schemes.Compact
        elif self.flux_mode == 'compact adaptive':
            self.flux_calculator = flux_schemes.CompactAdaptive
        elif self.flux_mode == 'slip':
            self.flux_calculator = flux_schemes.SLIP
        elif self.flux_mode == 'minmod':
            self.flux_calculator = flux_schemes.MinmodLimited
        else:
            raise RuntimeError
        self.flux_calculator = self.flux_calculator(self.N,
                                                    self.dz,
                                                    self.V,
                                                    self.z,
                                                    self.u)
        self.flux_calculator.update_ion_parameters(self.equlibrator)

    def set_reference_frame(self, frame=None, edge='right'):
        """Set the frame of reference.

        Frame should be an ion. Edge should be right or left.
        """
        self.frame = self.flux_calculator.frame = frame
        self.edge = self.flux_calculator.edge = edge

    def _decompose_state(self, state):
        """Decompose the state into X and concentrations."""
        x = state[:self.N]
        concentrations = state[self.N:].reshape(self.concentrations.shape)
        return (x, concentrations)

    def _compose_state(self, x, concentrations):
        """Compose X and concentrations into a state."""
        x = x.flatten()
        concentrations = concentrations.flatten()
        state = np.concatenate((x, concentrations))
        return state

    def _write_solution(self, t, state, full=True):
        """Write the current state to solutions."""
        (x, concentrations) = self._decompose_state(state)
        pH = self.equlibrator.pH
        ionic_strength = self.equlibrator.ionic_strength
        current_electrolyte = \
            Electrolyte(nodes=x, ions=self.ions,
                        concentrations=concentrations,
                        pH=pH, ionic_strength=ionic_strength,
                        voltage=self.V, current_density=self.flux_calculator.j,
                        )
        self.electromigration.add_electrolyte(t, current_electrolyte, full)

    def solve(self, tmax, dt=1, method='dopri5'):
        """Solve for a series of time points using an ODE solver."""
        self.solver = solver = integrate.ode(self._objective)

        solver.set_integrator(method, atol=self.atol, rtol=self.rtol)

        if solver._integrator.supports_solout:
            solver.set_solout(self._solout)
        else:
            warnings.warn("Solver doesn't support solout.")

        solver.set_initial_value(self._compose_state(self.x,
                                                     self.concentrations
                                                     ))

        while solver.successful() and solver.t < tmax:
            tnew = solver.t + dt
            if tnew > tmax:
                tnew = tmax
            solver.integrate(tnew)
            self._write_solution(solver.t, solver.y, full=False)
        (self.x, self.concentrations) = self._decompose_state(solver.y)

        if not solver.successful():
            print 'solver failed at time', solver.t

    def _solout(self, t, state):
        """Perform actions when a successful solution step is found."""
        (self.x, self.concentrations) = self._decompose_state(state)
        self.equlibrator.equilibrate(self.concentrations)
        self.flux_calculator.update_ion_parameters(self.equlibrator)
        self._write_solution(t, state, full=True)

    def _objective(self, t, state):
        """The objective function of the solver."""
        self.t = t
        (self.x, self.concentrations) = self._decompose_state(state)
        dcdt = self.flux_calculator.dcdt(self.x, self.concentrations)
        dxdt = self.flux_calculator.node_flux
        dstatedt = self._compose_state(dxdt, dcdt)
        return dstatedt
